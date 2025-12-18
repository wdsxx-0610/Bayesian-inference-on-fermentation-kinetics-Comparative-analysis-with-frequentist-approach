# ==============================================================================
# Task: Extract priors from Co-fermentation data (single file)
# ==============================================================================

graphics.off()
rm(list=ls())

# --- 1. Load packages ---
library(rjags)
library(coda)
library(tidyverse)
library(ggmcmc)

# --- 2. Path configuration ---
# Modify to your ribocose_cleaned.csv path
data_file_path <- "../data/ribocose_cleaned.csv"
save_dir <- "./priors_extraction_results"

if (!dir.exists(save_dir)) dir.create(save_dir)

# --- [Critical] Column name mapping table ---
# Format: "Standard metabolite name" = "Column name in CSV"
# According to your description: lactate, acetate, ethanol, PDO
col_map <- list(
  "Lactate" = "lactate",
  "Acetate" = "acetate",
  "Ethanol" = "ethanol",
  "PDO"     = "PDO"
)

# --- 3. Read main data file ---
if (!file.exists(data_file_path)) {
  stop(paste("Data file does not exist:", data_file_path))
}
main_df <- read.csv(data_file_path)
cat("Successfully read data file, contains columns:", colnames(main_df), "\n")

# --- 4. Define extraction function ---
extract_prior_from_column <- function(metabolite_name, col_name) {
  
  cat(paste0("\n==================================================\n"))
  cat(paste0("Processing: ", metabolite_name, " (column name: ", col_name, ")\n"))
  
  # 4.1 Data preparation
  # Check if column exists
  if (!col_name %in% colnames(main_df)) {
    warning(paste("Column name", col_name, "does not exist in data, skipping."))
    return(NULL)
  }
  
  # Extract Time (column 1) and Target Data
  # Assumes column 1 is time 'x' or 'time'
  xData <- main_df[, 1] 
  yData <- main_df[[col_name]]
  
  # Remove NA
  valid_idx <- !is.na(yData)
  xData <- xData[valid_idx]
  yData <- yData[valid_idx]
  
  cat(paste0("  -> Number of data points: ", length(yData), "\n"))
  
  # 4.2 Initial value guess (helps convergence)
  sdY <- sd(yData)
  guess_Amax <- max(yData)
  if(guess_Amax == 0) guess_Amax <- 0.1 # Prevent all-zero data error
  
  initsList <- list(
    Amax = guess_Amax,
    a = 0.1, 
    b = 0.001, 
    c = 10,    
    sigma = ifelse(sdY==0, 0.1, sdY) # Prevent sd=0
  )
  
  # 4.3 JAGS Model (Weak Priors)
  # Use exactly the same MHSF model structure, let data determine parameters
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
  
  # 4.4 Run JAGS
  dataList <- list(x = xData, y = yData, Ntotal = length(yData), init_Amax = guess_Amax)
  
  jagsModel <- jags.model("Temp_Model_Coferm.txt", data = dataList, inits = initsList, 
                          n.chains = 3, n.adapt = 1000, quiet = TRUE)
  update(jagsModel, n.iter = 2000)
  codaSamples <- coda.samples(jagsModel, variable.names = c("Amax", "a", "b", "c", "sigma"), 
                              n.iter = 5000, thin = 2)
  
  # 4.5 Diagnostics
  gelman_out <- gelman.diag(codaSamples, multivariate = FALSE)
  psrf <- gelman_out$psrf[,1]
  
  if(any(psrf > 1.1)) {
    cat("  ⚠️ Warning: Poor convergence (Gelman > 1.1)\n")
  } else {
    cat("  ✅ Good convergence\n")
  }
  
  # 4.6 Extract results (Mean & SD)
  sum_stats <- summary(codaSamples)$statistics
  prior_params <- data.frame(
    Metabolite = metabolite_name,
    Parameter = rownames(sum_stats),
    Mean = sum_stats[, "Mean"],
    SD = sum_stats[, "SD"]
  )
  
  # Save Traceplot
  ggs_object <- ggs(codaSamples)
  ggsave(file.path(save_dir, paste0(metabolite_name, "_Trace.png")), 
         ggs_traceplot(ggs_object), width=8, height=6)
  
  return(prior_params)
}

# --- 5. Batch execution ---
results_list <- list()

for (met_name in names(col_map)) {
  col_name <- col_map[[met_name]]
  
  # Use tryCatch to prevent error from one metabolite interrupting entire process
  try({
    res <- extract_prior_from_column(met_name, col_name)
    if (!is.null(res)) {
      results_list[[met_name]] <- res
    }
  })
}

# --- 6. Summary output ---
final_priors <- do.call(rbind, results_list)
rownames(final_priors) <- NULL

print("--------------------------------------------------")
print("Extraction results preview:")
print(head(final_priors))
print("--------------------------------------------------")

# Save CSV
output_csv <- file.path(save_dir, "Extracted_Priors_from_Cofermentation.csv")
write.csv(final_priors, output_csv, row.names = FALSE)

cat(paste0("\nProcessing complete! Results saved to: ", output_csv, "\n"))
