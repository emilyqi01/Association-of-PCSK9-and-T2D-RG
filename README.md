# Investigating the effects of coding \textit{PCSK9} gene variants on blood glucose levels and susceptibility to type 2 diabetes

This repository contains scripts and data preprocessing workflows for analyzing the association between PCSK9 variants and Type 2 Diabetes (T2D) risk, including random glucose levels (RG) and lipid traits. The analysis is conducted using various statistical methods, including single-variant and combined-variant approaches.

## Table of Contents

- [Project Structure](#project-structure)
- [Setup and Installation](#setup-and-installation)
- [Data Preprocessing](#data-preprocessing)
- [Single Variant Analysis](#single-variant-analysis)
- [Combined Variants Analysis](#combined-variants-analysis)
- [Contributing](#contributing)
- [License](#license)

## Project Structure

### 1. `data_preprocess/`

This folder contains scripts for preprocessing the raw genetic and phenotypic data.

- `deal_raw_data.sh`: Shell script to clean and preprocess the raw data files.
- `lipid_data_transformation_and_visualization.R`: R script for transforming lipid phenotype data and visualizing distributions.
- `process_and_visualize_genetic_data.R`: R script to process genetic data and generate visual summaries.

### 2. `Single_variant_analysis/`

This folder includes scripts for performing single-variant analyses, such as linear regression and logistic regression, to identify significant SNPs.

#### Subdirectory: `linear_regression_classification/`

- `plink_code.sh`: Shell script to run PLINK for linear regression analysis.
- `significant_variants.R`: R script to identify significant variants using logistic regression and classification methods.

#### Subdirectory: `regression_classification_penalty/`

- `encode_data.py`: Python script to encode data for machine learning models.
- `rg_tune_lasso.py`: Python script for tuning Lasso regression models for random glucose levels.
- `sgd_lasso_t2d.py`: Python script to apply stochastic gradient descent (SGD) with Lasso regularization for T2D analysis.

### 3. `Combined_variants_analysis/`

This folder contains scripts for combined variant analysis using burden tests, SKAT, and hierarchical clustering.

- `SAIGE_scripts.sh`: Shell script to run SAIGE for group-based analysis.
- `analysis_combined_pvalues.R`: R script to analyze combined p-values from multiple tests.
- `hierarchical_clustering.R`: R script for hierarchical clustering of genetic variants based on Z-scores.

## Setup and Installation

To run the scripts in this repository, you will need to have R, Python, and PLINK installed, along with the necessary libraries:

### R Libraries

- `ggplot2`, `dplyr`, `tidyr`, `gridExtra`, `MASS`

### Python Libraries

- `pandas`, `scikit-learn`, `numpy`

Install the libraries using `install.packages()` in R or `pip install` for Python.

## Data Preprocessing

Before running the analysis scripts, ensure that the raw data is processed correctly. Use the scripts in the `data_preprocess/` directory to clean and transform the data as needed.

\`\`\`bash
# Example command to preprocess data
bash data_preprocess/deal_raw_data.sh
\`\`\`

## Single Variant Analysis

To conduct single variant analysis, navigate to the `Single_variant_analysis/linear_regression_classification/` directory and run the scripts provided.

\`\`\`bash
# Run PLINK for single variant analysis
bash Single_variant_analysis/linear_regression_classification/plink_code.sh
\`\`\`

## Combined Variants Analysis

Combined variant analysis includes the use of statistical tests like Burden, SKAT, and SKAT-O. Navigate to the `Combined_variants_analysis/` directory and run the respective scripts.

\`\`\`bash
# Run SAIGE for combined variants analysis
bash Combined_variants_analysis/SAIGE_scripts.sh
\`\`\`

## Contributing

Contributions are welcome! Please fork this repository and submit a pull request for any proposed changes or additions.

## License
