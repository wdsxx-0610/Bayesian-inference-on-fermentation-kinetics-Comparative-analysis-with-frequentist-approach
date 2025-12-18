# ==============================================================================
# Task: Extract priors from Test A data (adapted for multi-file structure)
# ==============================================================================
install.packages("ggmcmc")
graphics.off()
rm(list=ls())

# --- 1. Load packages ---
library(rjags)
library(coda)
library(tidyverse)
library(ggmcmc) 

# --- 2. Path configuration ---
# Your data folder path
base_dir <- "./biohydrogendata"
save_dir <- "./priors_extraction_results"

if (!dir.exists(save_dir)) dir.create(save_dir)

# --- [Critical] File name mapping table ---
# Please ensure the file names here match the actual file names in your folder
# Format: "Metabolite name" = "filename.csv"
file_map <- list(
  "Acetate"  = "acetate.csv",
  "Hydrogen" = "hydrogen.csv",
  "Biomass"  = "biomass.csv",
  "Lactate"  = "lactate.csv",
  "Glucose"  = "glucose.csv" 
)

# Test A data column number (according to your description, time is x(1), then TestB(2), TestA is (3))
# If different files have different column numbers, please adjust individually, here assumes all are column 3
target_col_idx <- 3 

# --- 3. Define extraction function ---
extract_prior_from_file <- function(metabolite_name, file_name) {
  
  full_path <- file.path(base_dir, file_name)
  
  # Check if file exists
  if (!file.exists(full_path)) {
    warning(paste("File does not exist:", full_path))
    return(NULL)
  }
  
  cat(paste0("\n==================================================\n"))
  cat(paste0("Processing: ", metabolite_name, " (reading ", file_name, ")\n"))
  
  # 3.1 Read data
  df <- read.csv(full_path)
  
  # Extract Time (column 1) and Test A Data (column 3)
  # Use na.omit to remove null values to prevent errors
  clean_data <- df[, c(1, target_col_idx)] %>% na.omit()
  xData <- clean_data[, 1]
  yData <- clean_data[, 2]
  
  cat(paste0("  -> Number of data points: ", length(yData), "\n"))
  
  # 3.2 Initial value guess
  sdY <- sd(yData)
  guess_Amax <- max(yData)
  
  initsList <- list(
    Amax = guess_Amax,
    a = 0.1, 
    b = 0.001, 
    c = 10,    
    sigma = sdY
  )
  
  # 3.3 JAGS Model (Weak Priors)
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      y[i] ~ dnorm(mu[i], sigma_tau)
      mu[i] <- Amax / (exp(-a * (x[i] - c)) + exp(b * (x[i] - c)))
    }
    
    # Weak Priors - let Test A data determine parameters
    Amax ~ dnorm(init_Amax, 1.0E-4) T(0,) 
    a    ~ dlnorm(0, 1.0)  
    b    ~ dlnorm(0, 1.0)
    c    ~ dnorm(20, 1.0E-4)
    
    sigma ~ dunif(0, 100)
    sigma_tau <- 1 / (sigma * sigma)
  }
  "
  writeLines(modelString, con = "Temp_Model.txt")
  
  # 3.4 Run JAGS
  dataList <- list(x = xData, y = yData, Ntotal = length(yData), init_Amax = guess_Amax)
  
  jagsModel <- jags.model("Temp_Model.txt", data = dataList, inits = initsList, 
                          n.chains = 3, n.adapt = 1000, quiet = TRUE)
  update(jagsModel, n.iter = 2000)
  codaSamples <- coda.samples(jagsModel, variable.names = c("Amax", "a", "b", "c", "sigma"), 
                              n.iter = 10000, thin = 2)
  
  # 3.5 Diagnostics and results
  gelman_out <- gelman.diag(codaSamples, multivariate = FALSE)
  psrf <- gelman_out$psrf[,1]
  
  if(any(psrf > 1.1)) {
    cat("  ⚠️ Warning: Poor convergence (Gelman > 1.1)\n")
  } else {
    cat("  ✅ Good convergence\n")
  }
  
  # Extract statistics
  sum_stats <- summary(codaSamples)$statistics
  prior_params <- data.frame(
    Metabolite = metabolite_name,
    Parameter = rownames(sum_stats),
    Mean = sum_stats[, "Mean"],
    SD = sum_stats[, "SD"]
  )
  
  # Plot and save
  ggs_object <- ggs(codaSamples)
  ggsave(file.path(save_dir, paste0(metabolite_name, "_Trace.png")), 
         ggs_traceplot(ggs_object), width=8, height=6)
  
  return(prior_params)
}

# --- 4. Batch execution ---
results_list <- list()

for (met_name in names(file_map)) {
  file_name <- file_map[[met_name]]
  # Use tryCatch to prevent failure from one file causing interruption
  try({
    res <- extract_prior_from_file(met_name, file_name)
    if (!is.null(res)) {
      results_list[[met_name]] <- res
    }
  })
}

# --- 5. Summary output ---
final_priors <- do.call(rbind, results_list)
rownames(final_priors) <- NULL

print(final_priors)

# Save final results
write.csv(final_priors, file.path(save_dir, "Extracted_Priors_from_TestA.csv"), row.names = FALSE)
cat(paste0("\nProcessing complete! Results saved to: ", save_dir, "\n"))
