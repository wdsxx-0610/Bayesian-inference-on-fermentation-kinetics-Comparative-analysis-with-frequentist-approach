# ==============================================================================
# Script: Convert Test A natural parameters to Test C Log-Normal prior parameters
# Output: Generate a new CSV file with mu_log and tau_log ready to fill into JAGS
# ==============================================================================

rm(list=ls())
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Configure paths
# ------------------------------------------------------------------------------
# Input file: previously extracted Test A raw results
input_file  <- "../priors/priors_extraction_results/Extracted_Priors_from_TestA.csv"

# Output file: file prepared for JAGS
output_file <- "../priors/priors_extraction_results/Ready_for_JAGS_Priors.csv"

# Set engineering parameters
inflation_factor <- 2.0  # SD inflation coefficient (relax 2 times)

# ------------------------------------------------------------------------------
# 2. Define conversion function
# ------------------------------------------------------------------------------
# Convert natural scale Mean/SD to log scale Mean/Precision
# m: natural mean (already scaled)
# s: natural standard deviation (already scaled and inflated)
calc_log_params <- function(m, s) {
  # Prevent m or s from being 0 causing calculation errors
  m <- max(m, 1e-6)
  s <- max(s, 1e-6)
  
  # Formula
  var_log <- log(1 + (s^2 / m^2))
  mu_log  <- log(m) - 0.5 * var_log
  tau_log <- 1 / var_log
  
  return(list(mu_log = mu_log, tau_log = tau_log))
}

# ------------------------------------------------------------------------------
# 3. Read and process data
# ------------------------------------------------------------------------------
cat(">>> Reading data from:", input_file, "\n")
df <- read.csv(input_file)

# Create empty columns to store results
df$Scale_Factor <- NA
df$Target_Mean_Nat <- NA  # Target natural mean (Test C)
df$Target_SD_Nat   <- NA  # Target natural standard deviation (Test C)
df$Log_Mean        <- NA  # dlnorm parameter 1
df$Log_Precision   <- NA  # dlnorm parameter 2 (Tau)

cat(">>> Processing parameters...\n")

# Iterate through each row for calculation
for (i in 1:nrow(df)) {
  
  row <- df[i, ]
  metabolite <- row$Metabolite
  param <- row$Parameter
  mean_val <- row$Mean
  sd_val <- row$SD
  
  # --- Step A: Determine scaling factor ---
  # Hydrogen (yield) -> 1.0
  # Acetate/Glucose/Biomass (concentration) -> 2.0
  if (metabolite == "Hydrogen") {
    s_factor <- 1.0
  } else {
    # For Amax (concentration parameter) double it, other kinetic parameters (a,b,c) not doubled
    if (param == "Amax") {
      s_factor <- 2.0
    } else {
      s_factor <- 1.0 # Kinetic constants usually don't change mean with substrate concentration
    }
  }
  
  # --- Step B: Calculate target natural parameters (Test C Natural Scale) ---
  # Mean: multiply by scaling factor
  target_mean <- mean_val * s_factor
  
  # Standard deviation: first multiply by scaling factor (if concentration), then by inflation coefficient (relax)
  target_sd <- sd_val * s_factor * inflation_factor
  
  # --- Step C: Logarithmic transformation ---
  log_res <- calc_log_params(target_mean, target_sd)
  
  # --- Step D: Fill back into table ---
  df$Scale_Factor[i]    <- s_factor
  df$Target_Mean_Nat[i] <- round(target_mean, 4)
  df$Target_SD_Nat[i]   <- round(target_sd, 4)
  df$Log_Mean[i]        <- round(log_res$mu_log, 4)
  df$Log_Precision[i]   <- round(log_res$tau_log, 4)
}

# ------------------------------------------------------------------------------
# 4. Save results
# ------------------------------------------------------------------------------
write.csv(df, output_file, row.names = FALSE)

# Print preview
print("----------------------------------------------------------------")
print("Conversion complete! Below is a preview of some results (Log_Mean and Log_Precision can be directly filled into dlnorm):")
print("----------------------------------------------------------------")
print(head(df[, c("Metabolite", "Parameter", "Target_Mean_Nat", "Log_Mean", "Log_Precision")]))

cat(paste0("\nFile saved to: ", output_file, "\n"))