# Bayesian-inference-in-fermentation-compared-with-frequentist-method

This study explores the application of Bayesian inference to fermentation kinetic modeling, focusing on two representative batch systems: glycerol-glucose co-fermentation for 1,3-propanediol (1,3-PDO) production and dark fermentation for biohydrogen production.

## Repository Structure

The repository is organized into the following directories for better maintainability:

### Directories

- **`data/`**: Contains all CSV data files used in the analysis
  - `CofermentationData.csv`: Co-fermentation experimental data
  - `biohydrogenData.csv`: Biohydrogen fermentation data
  - `ribocose_cleaned.csv`: Cleaned ribose data for prior extraction

- **`priors/`**: Scripts for extracting and processing prior distributions
  - `Get_Priors_from_TestA.R`: Extract priors from Test A data (multi-file structure)
  - `Get_priors_form_ribo.R`: Extract priors from co-fermentation data

- **`cofermentation/`**: Scripts for co-fermentation analysis
  - `Jags-cofermentation1211.R`: Main JAGS analysis script for co-fermentation
  - `Jags-cofermentation-source1211.r`: JAGS model source for co-fermentation
  - `F_cofermentation_regression_model.r`: Frequentist regression model
  - `FvsBplot_cofermentation_8.5.R`: Frequentist vs Bayesian plotting

- **`biohydrogen/`**: Scripts for biohydrogen fermentation analysis
  - `Jags-biohydrogen_1211.R`: Main JAGS analysis script for biohydrogen
  - `Jags-biohydrogen-source_1211.r`: JAGS model source for biohydrogen
  - `biohydrogen_regression_model.r`: Frequentist regression model

- **`prediction/`**: Scripts for prediction and comparison
  - `HalfPrediction_fit_1213.R`: Half-data prediction fitting (Frequentist vs Bayesian)
  - `HalfPrediction_plot_1212.R`: Visualization of prediction results

- **`utils/`**: Utility functions and helper scripts
  - `DBDA2E-utilities.R`: Utility functions from "Doing Bayesian Data Analysis" book
  - `normtolog.R`: Convert normal parameters to log-normal priors

- **`biohydrogendata/`**: Additional biohydrogen data files (acetate, biomass, glucose, hydrogen, lactate)

## Usage

### Running Analysis Scripts

All scripts use relative paths from their respective directories. To run a script:

1. Set your working directory to the script's directory
2. Run the script (e.g., `source("Jags-cofermentation1211.R")`)

Results will be saved to `./results/` subdirectories within each analysis folder.

### Dependencies

The scripts require the following R packages:
- `rjags`: JAGS interface for R
- `coda`: Output analysis for MCMC
- `tidyverse`: Data manipulation and visualization
- `ggplot2`: Advanced plotting
- `minpack.lm`: Nonlinear least squares fitting
- `Metrics`: Model evaluation metrics
- `investr`: Prediction intervals
- `ggmcmc`: MCMC diagnostics
- `cowplot`: Plot arrangements

Install JAGS separately from: http://mcmc-jags.sourceforge.net/

### Workflow

1. **Prior Extraction** (optional): Extract priors from previous datasets
   - Run scripts in `priors/` directory
   
2. **Cofermentation Analysis**: 
   - Run scripts in `cofermentation/` directory for Bayesian and Frequentist analysis
   
3. **Biohydrogen Analysis**: 
   - Run scripts in `biohydrogen/` directory for Bayesian and Frequentist analysis
   
4. **Prediction Analysis**: 
   - Run scripts in `prediction/` directory to compare methods using half-data

## Notes

- All paths in scripts are now relative to support independent execution
- Results are saved in `./results/` subdirectories
- Scripts create necessary output directories automatically
- All comments have been translated to English for broader accessibility 
