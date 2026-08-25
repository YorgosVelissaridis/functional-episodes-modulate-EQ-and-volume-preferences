# Functional Episodes Modulate Equalization and Volume Preferences

This repository contains the data and R pipeline for the manuscript: **"Functional Episodes modulate equalization and volume preferences in simulated automotive music listening"**.

It features a Bayesian hierarchical modeling approach (using `brms` and Stan) to estimate how internal listening goals, car cabin noise, and the music being reproduced shape listener audio preferences, evaluated via Leave-One-Covariate-Out (LOCO) cross-validation.

## 📂 Repository Structure

* `data/`: Contains the experimental dataset (pre-processed: `data_en_screened.csv`).

* `R/`: Contains the complete analytic pipeline.

  * `00_config.R`: Global variables, model formulas, paths, and MCMC execution flags.

  * `01_main_models.R`: Fits the primary continuous (Volume) and categorical (EQ) models, extracts R², computes odds ratios, and renders the main multi-panel coefficient distribution plot.

  * `02_loco_analysis.R`: Computes exact out-of-sample predictive accuracy drops (ΔELPD) by refitting 20 reduced models, applying Pareto-smoothed importance sampling and Bayesian bootstrapping.

  * `03_supplement_sensitivity.R`: Computes individual-level episode sensitivity metrics (volume SD and EQ perplexity) and generates the behavioral distribution histograms.

  * `04_supplement_confounds.R`: Evaluates auxiliary confound models (mood, hearing, immersion, task difficulty/duration) via LOO cross-validation and models the volume-EQ coupling.

  * `utils.R`: Shared `ggplot2` and `ggdist` visualization wrappers.

* `results/`: Contains the final rendered tables (`.csv`) and plots (`.pdf`).

  * `cache/`: A local directory (ignored by git) where the computationally heavy `.rds` MCMC models and `loo` objects are saved to prevent redundant compilation.

## ⚙️ Prerequisites

To guarantee exact computational reproducibility, this project uses `renv` to manage package dependencies (including `brms`, `loo`, `tidyverse`, and `cmdstanr`).

Note: Because this project uses `cmdstanr` for highly efficient Bayesian sampling, you must have a working C++ toolchain and CmdStan installed on your system prior to running the models.

## 🚀 How to Run

1. **Clone the repository** and set your working directory to the project root.

2. **Restore the environment**: Open the project in R or RStudio. `renv` will automatically bootstrap itself. Run `renv::restore()` in the console to download and install the exact package versions used to generate the manuscript results. If you want to run the python files, you can create a virtual environment and use `requirements.txt` to install the required libraries.

3. *(Optional)* **Review configurations** in `R/00_config.R`.

4. **Run the files**: Execute any of the `R` or `python` files.

*Note: Bayesian categorical MCMC sampling and exhaustive LOO calculations are highly computationally intensive. The initial run of the full pipeline may take hours depending on your hardware, but all subsequent runs will load instantly from the local `results/cache/`.*