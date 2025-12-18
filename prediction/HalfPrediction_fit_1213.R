# ==============================================================================
# Complete nonlinear regression comparison code: Frequentist vs Bayesian (User Specified Priors)
# Modified: 2024-12-13
# ==============================================================================

# --- 1. Load required packages ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rjags, coda, Metrics, investr)

# --- 2. Set paths (modify according to your actual situation) ---
data_path <- "../data/CofermentationData.csv" 
save_dir  <- "./results/Half_prediction_Results"
if (!dir.exists(save_dir)) dir.create(save_dir)

# --- 3. Define function to get JAGS model for specific metabolite (using user-provided Log-Normal priors) ---
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

# --- 4. Core fitting function ---
fit_metabolite <- function(df, metabolite_col) {
  
  cat(sprintf("\n=========================================\n"))
  cat(sprintf("Processing Metabolite: %s", metabolite_col))
  cat(sprintf("\n=========================================\n"))
  
  # Data preparation
  if (!metabolite_col %in% colnames(df)) {
    stop(paste("Column", metabolite_col, "not found in dataframe."))
  }
  
  df_sub <- df %>% 
    select(time, !!sym(metabolite_col)) %>% 
    rename(value = !!sym(metabolite_col)) %>%
    filter(!is.na(value)) # Remove null values
  
  n <- nrow(df_sub)
  if(n == 0) stop("No data for this metabolite.")
  
  # 50% for training, 50% for testing
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
  
  # Prediction range
  full_time <- seq(0, 70, length.out = 200)
  
  if (!is.null(nonlinear_model)) {
    # Prediction interval
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
    
    # Test set RMSE
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
  
  # Initial values: can give a reasonable starting exp value based on log-normal mean
  # For simplicity, give a relatively general initial value, or use training data
  inits <- function() list(
    logAmax = log(max(train_df$value)),
    logc    = log(mean(train_df$time)), 
    logtau  = 0
  )
  
  model_string <- get_model_string(metabolite_col)
  writeLines(model_string, "temp_model.jags")
  
  # Run JAGS
  # Increase n.adapt and n.iter to ensure convergence, because strong priors are used
  model <- jags.model("temp_model.jags", data = data_jags, inits = inits, n.chains = 3, n.adapt = 2000, quiet = FALSE)
  update(model, 2000)
  # Monitor Amax, a, b, c, tau (physical parameters)
  samples <- coda.samples(model, variable.names = c("Amax", "a", "b", "c", "tau"), n.iter = 10000)
  post <- as.data.frame(do.call(rbind, samples))
  
  # Bayesian prediction interval
  pred_matrix <- sapply(full_time, function(tx) {
    # Vectorized calculation for acceleration
    mu <- post$Amax / (exp(-post$a * (tx - post$c)) + exp(post$b * (tx - post$c)))
    sigma <- 1 / sqrt(post$tau)
    rnorm(length(mu), mu, sigma)
  })
  
  pred_bayes  <- apply(pred_matrix, 2, mean)
  lower_bayes <- apply(pred_matrix, 2, quantile, 0.025)
  upper_bayes <- apply(pred_matrix, 2, quantile, 0.975)
  
  # Test set RMSE
  test_pred_bayes <- sapply(test_df$time, function(tx) {
    mean(post$Amax / (exp(-post$a * (tx - post$c)) + exp(post$b * (tx - post$c))))
  })
  
  rmse_bayes <- rmse(test_df$value, test_pred_bayes)
  nsme_bayes <- 1 - sum((test_pred_bayes - test_df$value)^2) / sum((test_df$value - mean(test_df$value))^2)
  
  # ---------------------------------------------------------
  # C. Save results
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

# --- 5. Main execution logic (modify here to run separately) ---

# --- 5.1 Read and clean data ---
if (file.exists(data_path)) {
  raw_df <- read.csv(data_path)
  
  # Clean column names: 
  # 1. Ensure x -> time
  # 2. Capitalize first letter of other columns (so lactate becomes Lactate)
  df <- raw_df %>%
    rename(time = time) %>%
    rename_with(~ str_to_title(.), .cols = -time) 
  
  # [Critical fix] If auto capitalization changed PDO to Pdo, force it back to PDO
  if ("Pdo" %in% colnames(df)) {
    df <- df %>% rename(PDO = Pdo)
  }
  
  print("Data loaded successfully from CofermentationData.csv")
  print("Columns found:")
  print(colnames(df)) # Print column names, confirm "PDO" is present
  print(head(df))
  
} else {
  stop(paste("File not found:", data_path))
}

# 5.2 Select which metabolite to run
# ==========================================================
# Modify the target_metabolite variable below to select which one to run
# Options: "Lactate", "Acetate", "Ethanol", "PDO"
# ==========================================================

target_metabolite <- "Lactate"   # <--- Modify here! For example, change to "Acetate"
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
