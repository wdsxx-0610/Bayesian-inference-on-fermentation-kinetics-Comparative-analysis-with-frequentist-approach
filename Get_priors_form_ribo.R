# ==============================================================================
# 任务：从 Co-fermentation 数据 (单一文件) 中提取先验
# ==============================================================================

graphics.off()
rm(list=ls())

# --- 1. 加载包 ---
library(rjags)
library(coda)
library(tidyverse)
library(ggmcmc)

# --- 2. 路径配置 ---
# 请修改为你的 ribocose_cleaned.csv 所在路径
data_file_path <- "/Users/wdsxx0610/Documents/R_directory/A_formation/cofermentation/ribocose_cleaned.csv"
save_dir <- "/Users/wdsxx0610/Documents/R_directory/A_formation/cofermentation/"

if (!dir.exists(save_dir)) dir.create(save_dir)

# --- [关键] 列名映射表 ---
# 格式: "标准代谢物名称" = "CSV中的列名"
# 根据你的描述：lactate, acetate, ethanol, PDO
col_map <- list(
  "Lactate" = "lactate",
  "Acetate" = "acetate",
  "Ethanol" = "ethanol",
  "PDO"     = "PDO"
)

# --- 3. 读取主数据文件 ---
if (!file.exists(data_file_path)) {
  stop(paste("数据文件不存在:", data_file_path))
}
main_df <- read.csv(data_file_path)
cat("成功读取数据文件，包含列:", colnames(main_df), "\n")

# --- 4. 定义提取函数 ---
extract_prior_from_column <- function(metabolite_name, col_name) {
  
  cat(paste0("\n==================================================\n"))
  cat(paste0("正在处理: ", metabolite_name, " (列名: ", col_name, ")\n"))
  
  # 4.1 数据准备
  # 检查列是否存在
  if (!col_name %in% colnames(main_df)) {
    warning(paste("列名", col_name, "在数据中不存在，跳过。"))
    return(NULL)
  }
  
  # 提取 Time (第1列) 和 Target Data
  # 假设第1列是时间 'x' 或 'time'
  xData <- main_df[, 1] 
  yData <- main_df[[col_name]]
  
  # 去除 NA
  valid_idx <- !is.na(yData)
  xData <- xData[valid_idx]
  yData <- yData[valid_idx]
  
  cat(paste0("  -> 数据点数量: ", length(yData), "\n"))
  
  # 4.2 初始值猜测 (帮助收敛)
  sdY <- sd(yData)
  guess_Amax <- max(yData)
  if(guess_Amax == 0) guess_Amax <- 0.1 # 防止全0数据报错
  
  initsList <- list(
    Amax = guess_Amax,
    a = 0.1, 
    b = 0.001, 
    c = 10,    
    sigma = ifelse(sdY==0, 0.1, sdY) # 防止 sd=0
  )
  
  # 4.3 JAGS 模型 (Weak Priors)
  # 使用完全相同的 MHSF 模型结构，让数据决定参数
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dnorm(mu[i], sigma_tau)
      mu[i] <- Amax / (exp(-a * (x[i] - c)) + exp(b * (x[i] - c)))
    }
    
    # Weak Priors
    Amax ~ dnorm(init_Amax, 1.0E-4) T(0,) 
    a    ~ dlnorm(0, 1.0)  
    b    ~ dlnorm(0, 1.0)
    c    ~ dnorm(20, 1.0E-4)
    
    sigma ~ dunif(0, 100)
    sigma_tau <- 1 / (sigma * sigma)
  }
  "
  writeLines(modelString, con = "Temp_Model_Coferm.txt")
  
  # 4.4 运行 JAGS
  dataList <- list(x = xData, y = yData, Ntotal = length(yData), init_Amax = guess_Amax)
  
  jagsModel <- jags.model("Temp_Model_Coferm.txt", data = dataList, inits = initsList, 
                          n.chains = 3, n.adapt = 1000, quiet = TRUE)
  update(jagsModel, n.iter = 2000)
  codaSamples <- coda.samples(jagsModel, variable.names = c("Amax", "a", "b", "c", "sigma"), 
                              n.iter = 5000, thin = 2)
  
  # 4.5 诊断
  gelman_out <- gelman.diag(codaSamples, multivariate = FALSE)
  psrf <- gelman_out$psrf[,1]
  
  if(any(psrf > 1.1)) {
    cat("  ⚠️ 警告: 收敛性不佳 (Gelman > 1.1)\n")
  } else {
    cat("  ✅ 收敛良好\n")
  }
  
  # 4.6 提取结果 (Mean & SD)
  sum_stats <- summary(codaSamples)$statistics
  prior_params <- data.frame(
    Metabolite = metabolite_name,
    Parameter = rownames(sum_stats),
    Mean = sum_stats[, "Mean"],
    SD = sum_stats[, "SD"]
  )
  
  # 保存 Traceplot
  ggs_object <- ggs(codaSamples)
  ggsave(file.path(save_dir, paste0(metabolite_name, "_Trace.png")), 
         ggs_traceplot(ggs_object), width=8, height=6)
  
  return(prior_params)
}

# --- 5. 批量执行 ---
results_list <- list()

for (met_name in names(col_map)) {
  col_name <- col_map[[met_name]]
  
  # 使用 tryCatch 防止某个代谢物报错中断整个流程
  try({
    res <- extract_prior_from_column(met_name, col_name)
    if (!is.null(res)) {
      results_list[[met_name]] <- res
    }
  })
}

# --- 6. 汇总输出 ---
final_priors <- do.call(rbind, results_list)
rownames(final_priors) <- NULL

print("--------------------------------------------------")
print("提取结果预览:")
print(head(final_priors))
print("--------------------------------------------------")

# 保存 CSV
output_csv <- file.path(save_dir, "Extracted_Priors_from_Cofermentation.csv")
write.csv(final_priors, output_csv, row.names = FALSE)

cat(paste0("\n处理完成！结果已保存至: ", output_csv, "\n"))
