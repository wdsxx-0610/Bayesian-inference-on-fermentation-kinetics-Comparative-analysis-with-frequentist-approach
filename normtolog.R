# ==============================================================================
# 脚本：将 Test A 的自然参数转换为 Test C 的 Log-Normal 先验参数
# 输出：生成一个新的 CSV 文件，包含直接可填入 JAGS 的 mu_log 和 tau_log
# ==============================================================================

rm(list=ls())
library(dplyr)

# ------------------------------------------------------------------------------
# 1. 配置路径
# ------------------------------------------------------------------------------
# 输入文件：之前提取的 Test A 原始结果
input_file  <- "/Users/wdsxx0610/Documents/R_directory/TestABdata/Priors_Extraction/Extracted_Priors_from_TestA.csv"

# 输出文件：准备好给 JAGS 用的文件
output_file <- "/Users/wdsxx0610/Documents/R_directory/TestABdata/Priors_Extraction/Ready_for_JAGS_Priors.csv"

# 设定工程参数
inflation_factor <- 2.0  # SD 膨胀系数 (放宽 2 倍)

# ------------------------------------------------------------------------------
# 2. 定义转换函数
# ------------------------------------------------------------------------------
# 将自然尺度的 Mean/SD 转换为 对数尺度的 Mean/Precision
# m: 自然均值 (已缩放)
# s: 自然标准差 (已缩放且膨胀)
calc_log_params <- function(m, s) {
  # 防止 m 或 s 为 0 导致计算错误
  m <- max(m, 1e-6)
  s <- max(s, 1e-6)
  
  # 公式
  var_log <- log(1 + (s^2 / m^2))
  mu_log  <- log(m) - 0.5 * var_log
  tau_log <- 1 / var_log
  
  return(list(mu_log = mu_log, tau_log = tau_log))
}

# ------------------------------------------------------------------------------
# 3. 读取并处理数据
# ------------------------------------------------------------------------------
cat(">>> Reading data from:", input_file, "\n")
df <- read.csv(input_file)

# 创建空的列来存储结果
df$Scale_Factor <- NA
df$Target_Mean_Nat <- NA  # 目标自然均值 (Test C)
df$Target_SD_Nat   <- NA  # 目标自然标准差 (Test C)
df$Log_Mean        <- NA  # dlnorm 的参数 1
df$Log_Precision   <- NA  # dlnorm 的参数 2 (Tau)

cat(">>> Processing parameters...\n")

# 遍历每一行进行计算
for (i in 1:nrow(df)) {
  
  row <- df[i, ]
  metabolite <- row$Metabolite
  param <- row$Parameter
  mean_val <- row$Mean
  sd_val <- row$SD
  
  # --- Step A: 确定缩放因子 ---
  # Hydrogen (产率) -> 1.0
  # Acetate/Glucose/Biomass (浓度) -> 2.0
  if (metabolite == "Hydrogen") {
    s_factor <- 1.0
  } else {
    # 对于 Amax (浓度参数) 进行翻倍，其他动力学参数(a,b,c) 不翻倍
    if (param == "Amax") {
      s_factor <- 2.0
    } else {
      s_factor <- 1.0 # 动力学常数通常不随底物浓度改变均值
    }
  }
  
  # --- Step B: 计算目标自然参数 (Test C Natural Scale) ---
  # 均值：乘以缩放因子
  target_mean <- mean_val * s_factor
  
  # 标准差：先乘以缩放因子(如果是浓度)，再乘以膨胀系数(放宽)
  target_sd <- sd_val * s_factor * inflation_factor
  
  # --- Step C: 对数转换 ---
  log_res <- calc_log_params(target_mean, target_sd)
  
  # --- Step D: 填回表格 ---
  df$Scale_Factor[i]    <- s_factor
  df$Target_Mean_Nat[i] <- round(target_mean, 4)
  df$Target_SD_Nat[i]   <- round(target_sd, 4)
  df$Log_Mean[i]        <- round(log_res$mu_log, 4)
  df$Log_Precision[i]   <- round(log_res$tau_log, 4)
}

# ------------------------------------------------------------------------------
# 4. 保存结果
# ------------------------------------------------------------------------------
write.csv(df, output_file, row.names = FALSE)

# 打印预览
print("----------------------------------------------------------------")
print("转换完成！以下是部分结果预览 (Log_Mean 和 Log_Precision 可直接填入 dlnorm):")
print("----------------------------------------------------------------")
print(head(df[, c("Metabolite", "Parameter", "Target_Mean_Nat", "Log_Mean", "Log_Precision")]))

cat(paste0("\n文件已保存至: ", output_file, "\n"))