# ==============================================================================
# 任务：从 Test A 数据中提取先验 (适配多文件结构)
# ==============================================================================
install.packages("ggmcmc")
graphics.off()
rm(list=ls())

# --- 1. 加载包 ---
library(rjags)
library(coda)
library(tidyverse)
library(ggmcmc) 

# --- 2. 路径配置 ---
# 你的数据所在的文件夹路径
base_dir <- "/Users/wdsxx0610/Documents/R_directory/TestABdata"
save_dir <- "/Users/wdsxx0610/Documents/R_directory/TestABdata/Priors_Extraction"

if (!dir.exists(save_dir)) dir.create(save_dir)

# --- [关键] 文件名映射表 ---
# 请确保这里的文件名与你文件夹里的实际文件名一致
# 格式: "代谢物名称" = "文件名.csv"
file_map <- list(
  "Acetate"  = "acetate.csv",
  "Hydrogen" = "hydrogen.csv",
  "Biomass"  = "biomass.csv",
  "Lactate"  = "lactate.csv",
  "Glucose"  = "glucose.csv" 
)

# Test A 数据所在的列号 (根据你的描述，时间是x(1)，然后是TestB(2)，TestA是(3))
# 如果不同文件列号不一样，请单独调整，这里假设都是第 3 列
target_col_idx <- 3 

# --- 3. 定义提取函数 ---
extract_prior_from_file <- function(metabolite_name, file_name) {
  
  full_path <- file.path(base_dir, file_name)
  
  # 检查文件是否存在
  if (!file.exists(full_path)) {
    warning(paste("文件不存在:", full_path))
    return(NULL)
  }
  
  cat(paste0("\n==================================================\n"))
  cat(paste0("正在处理: ", metabolite_name, " (读取 ", file_name, ")\n"))
  
  # 3.1 读取数据
  df <- read.csv(full_path)
  
  # 提取 Time (第1列) 和 Test A Data (第3列)
  # 使用 na.omit 去除空值，防止报错
  clean_data <- df[, c(1, target_col_idx)] %>% na.omit()
  xData <- clean_data[, 1]
  yData <- clean_data[, 2]
  
  cat(paste0("  -> 数据点数量: ", length(yData), "\n"))
  
  # 3.2 初始值猜测
  sdY <- sd(yData)
  guess_Amax <- max(yData)
  
  initsList <- list(
    Amax = guess_Amax,
    a = 0.1, 
    b = 0.001, 
    c = 10,    
    sigma = sdY
  )
  
  # 3.3 JAGS 模型 (Weak Priors)
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dnorm(mu[i], sigma_tau)
      mu[i] <- Amax / (exp(-a * (x[i] - c)) + exp(b * (x[i] - c)))
    }
    
    # Weak Priors - 让 Test A 数据决定参数
    Amax ~ dnorm(init_Amax, 1.0E-4) T(0,) 
    a    ~ dlnorm(0, 1.0)  
    b    ~ dlnorm(0, 1.0)
    c    ~ dnorm(20, 1.0E-4)
    
    sigma ~ dunif(0, 100)
    sigma_tau <- 1 / (sigma * sigma)
  }
  "
  writeLines(modelString, con = "Temp_Model.txt")
  
  # 3.4 运行 JAGS
  dataList <- list(x = xData, y = yData, Ntotal = length(yData), init_Amax = guess_Amax)
  
  jagsModel <- jags.model("Temp_Model.txt", data = dataList, inits = initsList, 
                          n.chains = 3, n.adapt = 1000, quiet = TRUE)
  update(jagsModel, n.iter = 2000)
  codaSamples <- coda.samples(jagsModel, variable.names = c("Amax", "a", "b", "c", "sigma"), 
                              n.iter = 10000, thin = 2)
  
  # 3.5 诊断与结果
  gelman_out <- gelman.diag(codaSamples, multivariate = FALSE)
  psrf <- gelman_out$psrf[,1]
  
  if(any(psrf > 1.1)) {
    cat("  ⚠️ 警告: 收敛性不佳 (Gelman > 1.1)\n")
  } else {
    cat("  ✅ 收敛良好\n")
  }
  
  # 提取统计量
  sum_stats <- summary(codaSamples)$statistics
  prior_params <- data.frame(
    Metabolite = metabolite_name,
    Parameter = rownames(sum_stats),
    Mean = sum_stats[, "Mean"],
    SD = sum_stats[, "SD"]
  )
  
  # 绘图保存
  ggs_object <- ggs(codaSamples)
  ggsave(file.path(save_dir, paste0(metabolite_name, "_Trace.png")), 
         ggs_traceplot(ggs_object), width=8, height=6)
  
  return(prior_params)
}

# --- 4. 批量执行 ---
results_list <- list()

for (met_name in names(file_map)) {
  file_name <- file_map[[met_name]]
  # 使用 tryCatch 防止某个文件读取失败导致中断
  try({
    res <- extract_prior_from_file(met_name, file_name)
    if (!is.null(res)) {
      results_list[[met_name]] <- res
    }
  })
}

# --- 5. 汇总输出 ---
final_priors <- do.call(rbind, results_list)
rownames(final_priors) <- NULL

print(final_priors)

# 保存最终结果
write.csv(final_priors, file.path(save_dir, "Extracted_Priors_from_TestA.csv"), row.names = FALSE)
cat(paste0("\n处理完成！结果已保存至: ", save_dir, "\n"))
