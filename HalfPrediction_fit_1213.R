# ==============================================================================
# 完整非线性回归对比代码：Frequentist vs Bayesian (User Specified Priors)
# Modified: 2024-12-13
# ==============================================================================

# --- 1. 加载必要的包 ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rjags, coda, Metrics, investr)

# --- 2. 设置路径 (请根据你的实际情况修改) ---
data_path <- "/Users/wdsxx0610/Documents/R_directory/A_formation/cofermentation/CofermentationData.csv" 
save_dir  <- "/Users/wdsxx0610/Documents/R_directory/A_formation/cofermentation/Half_prediction_Results"
if (!dir.exists(save_dir)) dir.create(save_dir)

# --- 3. 定义获取特定代谢物 JAGS 模型的函数 (使用用户提供的 Log-Normal 先验) ---
get_model_string <- function(metabolite_name) {
  
  if (metabolite_name == "Acetate") {
    # Prior (1): Acetate
    # inits <- c(4.0184, -1.7612, -5.7387, 3.2522, -0.012)
    # inits_precision <-c(46.6694/4, 27.2778/4, 1.4055/4, 120.2169/4, 1.4442)
    return("
    model {
      for (i in 1:n) {
        y[i] ~ dnorm(mu[i], tau)
        mu[i] <- Amax / (exp(-a * (x[i] - c)) + exp(b * (x[i] - c)))
      }
      
      # Priors on Log Scale (Log-Normal formulation)
      logAmax ~ dnorm(4.0184, 46.6694/4)
      loga    ~ dnorm(-1.7612, 27.2778/4)
      logb    ~ dnorm(-5.7387, 1.4055/4)
      logc    ~ dnorm(3.2522, 120.2169/4)
      logtau  ~ dnorm(-0.012, 1.4442)
      
      # Transform back to physical parameters
      Amax <- exp(logAmax)
      a    <- exp(loga)
      b    <- exp(logb)
      c    <- exp(logc)
      tau  <- exp(logtau)
      
      sigma <- 1 / sqrt(tau)
    }
    ")
    
  } else if (metabolite_name == "Ethanol") {
    # Prior (2): Ethanol
    # inits <- c(0.4095, -1.0635, -5.5425, 1.8806, -2.3111)
    # inits_precision <- c(0.3211, 0.8844, 0.2121, 0.6064, 0.3676)
    return("
    model {
      for (i in 1:n) {
        y[i] ~ dnorm(mu[i], tau)
        mu[i] <- Amax / (exp(-a * (x[i] - c)) + exp(b * (x[i] - c)))
      }
      
      # Priors on Log Scale
      logAmax ~ dnorm(0.4095, 0.3211)
      loga    ~ dnorm(-1.0635, 0.8844)
      logb    ~ dnorm(-5.5425, 0.2121)
      logc    ~ dnorm(1.8806, 0.6064)
      logtau  ~ dnorm(-2.3111, 0.3676)
      
      # Transform
      Amax <- exp(logAmax)
      a    <- exp(loga)
      b    <- exp(logb)
      c    <- exp(logc)
      tau  <- exp(logtau)
      
      sigma <- 1 / sqrt(tau)
    }
    ")
    
  } else if (metabolite_name == "Lactate") {
    # Prior (3): Lactate
    # inits <- c(2.7348, -1.6316, -1.7737, 3.0193, -0.6804)
    # inits_precision <-c(52.1589, 4.2834, 4.4571, 38.8296, 3.1878)
    return("
    model {
      for (i in 1:n) {
        y[i] ~ dnorm(mu[i], tau)
        mu[i] <- Amax / (exp(-a * (x[i] - c)) + exp(b * (x[i] - c)))
      }
      
      # Priors on Log Scale
      logAmax ~ dnorm(2.7348, 52.1589)
      loga    ~ dnorm(-1.6316, 4.2834)
      logb    ~ dnorm(-1.7737, 4.4571)
      logc    ~ dnorm(3.0193, 38.8296)
      logtau  ~ dnorm(-0.6804, 3.1878)
      
      # Transform
      Amax <- exp(logAmax)
      a    <- exp(loga)
      b    <- exp(logb)
      c    <- exp(logc)
      tau  <- exp(logtau)
      
      sigma <- 1 / sqrt(tau)
    }
    ")
    
  } else if (metabolite_name == "PDO") {
    # Prior (4): PDO
    # inits <- c(3.9939, -1.2502, -6.0002, 3.4651, 0.1206)
    # inits_precision <-c(25.3153, 12.7001, 0.7765, 202.8906, 0.9617)
    return("
    model {
      for (i in 1:n) {
        y[i] ~ dnorm(mu[i], tau)
        mu[i] <- Amax / (exp(-a * (x[i] - c)) + exp(b * (x[i] - c)))
      }
      
      # Priors on Log Scale
      logAmax ~ dnorm(3.9939, 25.3153)
      loga    ~ dnorm(-1.2502, 12.7001)
      logb    ~ dnorm(-6.0002, 0.7765)
      logc    ~ dnorm(3.4651, 202.8906)
      logtau  ~ dnorm(0.1206, 0.9617)
      
      # Transform
      Amax <- exp(logAmax)
      a    <- exp(loga)
      b    <- exp(logb)
      c    <- exp(logc)
      tau  <- exp(logtau)
      
      sigma <- 1 / sqrt(tau)
    }
    ")
    
  } else {
    stop(paste("Unknown metabolite:", metabolite_name))
  }
}

# --- 4. 核心拟合函数 ---
fit_metabolite <- function(df, metabolite_col) {
  
  cat(sprintf("\n=========================================\n"))
  cat(sprintf("Processing Metabolite: %s", metabolite_col))
  cat(sprintf("\n=========================================\n"))
  
  # 数据准备
  if (!metabolite_col %in% colnames(df)) {
    stop(paste("Column", metabolite_col, "not found in dataframe."))
  }
  
  df_sub <- df %>% 
    select(time, !!sym(metabolite_col)) %>% 
    rename(value = !!sym(metabolite_col)) %>%
    filter(!is.na(value)) # 去除空值
  
  n <- nrow(df_sub)
  if(n == 0) stop("No data for this metabolite.")
  
  # 50% 用于训练，50% 用于测试
  train_df <- df_sub[1:(n/2), ]
  test_df  <- df_sub[(n/2 + 1):n, ]
  
  # ---------------------------------------------------------
  # A. Frequentist Method (NLS)
  # ---------------------------------------------------------
  cat("  -> [1/2] Running Frequentist NLS...\n")
  
  nonlinear_model <- tryCatch({
    nls(
      value ~ Amax / (exp(-a * (time - c)) + exp(b * (time - c))),
      data = train_df,
      start = list(Amax = max(train_df$value), a = 0.2, b = 0.01, c = mean(train_df$time)),
      control = nls.control(maxiter = 500, warnOnly = TRUE)
    )
  }, error = function(e) {
    message("    NLS Fit failed: ", e$message)
    return(NULL)
  })
  
  # 预测范围
  full_time <- seq(0, 70, length.out = 200)
  
  if (!is.null(nonlinear_model)) {
    # 预测区间
    pred_fit_obj <- tryCatch({
      investr::predFit(
        nonlinear_model, 
        newdata = data.frame(time = full_time), 
        interval = "prediction", 
        level = 0.95
      )
    }, error = function(e) {
      # Fallback if investr fails
      matrix(NA, nrow=length(full_time), ncol=3, dimnames=list(NULL, c("fit","lwr","upr")))
    })
    
    train_pred_nls <- pred_fit_obj[, "fit"]
    ci_lower_nls   <- pred_fit_obj[, "lwr"]
    ci_upper_nls   <- pred_fit_obj[, "upr"]
    
    # 测试集RMSE
    test_pred_nls <- predict(nonlinear_model, newdata = test_df)
    rmse_nls <- rmse(test_df$value, test_pred_nls)
    nsme_nls <- 1 - sum((test_pred_nls - test_df$value)^2) / sum((test_df$value - mean(test_df$value))^2)
  } else {
    train_pred_nls <- rep(NA, length(full_time))
    ci_lower_nls <- rep(NA, length(full_time))
    ci_upper_nls <- rep(NA, length(full_time))
    rmse_nls <- NA; nsme_nls <- NA
  }
  
  # ---------------------------------------------------------
  # B. Bayesian Method (JAGS)
  # ---------------------------------------------------------
  cat("  -> [2/2] Running Bayesian JAGS with Specific Priors...\n")
  
  data_jags <- list(x = train_df$time, y = train_df$value, n = nrow(train_df))
  
  # 初始值: 可以根据 log-normal 的均值给一个合理的起始 exp 值
  # 这里为了简便，给一个相对通用的初始值，或者利用训练数据
  inits <- function() list(
    logAmax = log(max(train_df$value)),
    logc    = log(mean(train_df$time)), 
    logtau  = 0
  )
  
  model_string <- get_model_string(metabolite_col)
  writeLines(model_string, "temp_model.jags")
  
  # 运行 JAGS
  # 增加 n.adapt 和 n.iter 以确保收敛，因为使用了较强的先验
  model <- jags.model("temp_model.jags", data = data_jags, inits = inits, n.chains = 3, n.adapt = 2000, quiet = FALSE)
  update(model, 2000)
  # 监控 Amax, a, b, c, tau (物理参数)
  samples <- coda.samples(model, variable.names = c("Amax", "a", "b", "c", "tau"), n.iter = 10000)
  post <- as.data.frame(do.call(rbind, samples))
  
  # 贝叶斯预测区间
  pred_matrix <- sapply(full_time, function(tx) {
    # 向量化计算以加速
    mu <- post$Amax / (exp(-post$a * (tx - post$c)) + exp(post$b * (tx - post$c)))
    sigma <- 1 / sqrt(post$tau)
    rnorm(length(mu), mu, sigma)
  })
  
  pred_bayes  <- apply(pred_matrix, 2, mean)
  lower_bayes <- apply(pred_matrix, 2, quantile, 0.025)
  upper_bayes <- apply(pred_matrix, 2, quantile, 0.975)
  
  # 测试集 RMSE
  test_pred_bayes <- sapply(test_df$time, function(tx) {
    mean(post$Amax / (exp(-post$a * (tx - post$c)) + exp(post$b * (tx - post$c))))
  })
  
  rmse_bayes <- rmse(test_df$value, test_pred_bayes)
  nsme_bayes <- 1 - sum((test_pred_bayes - test_df$value)^2) / sum((test_df$value - mean(test_df$value))^2)
  
  # ---------------------------------------------------------
  # C. 保存结果
  # ---------------------------------------------------------
  save_df <- tibble(
    time = full_time,
    pred_nls = train_pred_nls,
    lower_nls = ci_lower_nls,
    upper_nls = ci_upper_nls,
    pred_bayes = pred_bayes,
    lower_bayes = lower_bayes,
    upper_bayes = upper_bayes
  )
  
  write_csv(save_df, file.path(save_dir, paste0(metabolite_col, "_predictions.csv")))
  
  cat(sprintf("  -> Done. RMSE NLS: %.4f | RMSE Bayes: %.4f\n", rmse_nls, rmse_bayes))
  
  return(list(
    Metabolite = metabolite_col,
    RMSE_Frequentist = rmse_nls, NSME_Frequentist = nsme_nls,
    RMSE_Bayesian = rmse_bayes, NSME_Bayesian = nsme_bayes
  ))
}

# --- 5. 主执行逻辑 (在此处修改以分别运行) ---

# --- 5.1 读取并清洗数据 ---
if (file.exists(data_path)) {
  raw_df <- read.csv(data_path)
  
  # 清洗列名: 
  # 1. 确保 x -> time
  # 2. 其他列首字母大写 (这样 lactate 会变成 Lactate)
  df <- raw_df %>%
    rename(time = time) %>%
    rename_with(~ str_to_title(.), .cols = -time) 
  
  # [关键修复] 如果自动大写把 PDO 变成了 Pdo，这里强制改回 PDO
  if ("Pdo" %in% colnames(df)) {
    df <- df %>% rename(PDO = Pdo)
  }
  
  print("Data loaded successfully from CofermentationData.csv")
  print("Columns found:")
  print(colnames(df)) # 打印列名，确认里面有 "PDO"
  print(head(df))
  
} else {
  stop(paste("找不到文件:", data_path))
}

# 5.2 选择要运行的代谢物
# ==========================================================
# 修改下方的 target_metabolite 变量来选择你要跑哪一个
# 可选项: "Lactate", "Acetate", "Ethanol", "PDO"
# ==========================================================

target_metabolite <- "Lactate"   # <--- 在这里修改! 例如改成 "Acetate"
target_metabolite <- "Acetate"  
target_metabolite <- "Ethanol"  
target_metabolite <- "PDO"  
# 5.3 运行逻辑
if (target_metabolite %in% colnames(df)) {
  
  # 运行单个
  result <- fit_metabolite(df, target_metabolite)
  
  # 显示结果
  print(as.data.frame(result))
  
  # 保存单个结果指标
  metrics_file <- file.path(save_dir, paste0(target_metabolite, "_metrics.csv"))
  write_csv(as.data.frame(result), metrics_file)
  cat(paste("\nMetrics saved to:", metrics_file, "\n"))
  
} else {
  stop(paste("Metabolite", target_metabolite, "not found in data columns (Check capitalization!)."))
}
