# Investigating the effects of coding *PCSK9* gene variants on blood glucose levels and susceptibility to type 2 diabetes

This repository contains scripts and data preprocessing workflows for analyzing the association between PCSK9 variants and Type 2 Diabetes (T2D) risk, including random glucose levels (RG) and lipid traits. The analysis is conducted using various statistical methods, including single-variant and combined-variant approaches.

## Table of Contents

- [Important Information Before Reproducing the Project](#important-information-before-reproducing-the-project)
- [Setup and Installation](#setup-and-installation)
- [Project Structure](#project-structure)
- [Data Preprocessing](#data-preprocessing)
- [Single Variant Analysis](#single-variant-analysis)
- [Combined Variants Analysis](#combined-variants-analysis)
- [Contributing](#contributing)
- [License](#license)

## Important Information Before Reproducing the Project

- **Data Access**: This project utilizes data from the UK Biobank, which is subject to copyright and cannot be redistributed. You will need to request access and download the necessary data directly from the [UK Biobank website](https://www.ukbiobank.ac.uk/).

- **Software Requirements**: The project relies on various software tools and is primarily executed via terminal commands. To successfully reproduce the results, you may need to install additional software. Please refer to the [Setup and Installation](#setup-and-installation) section for detailed instructions on the required software and installation steps.

- **Modifying Directory Paths**: The scripts in this project contain hardcoded directory paths. Before running the scripts, please update these paths to match your local environment.

## Setup and Installation

To run the analysis, you will need to install the SAIGE and PLINK software. It is recommended to use a virtual machine or a server to perform the analysis due to the computational requirements. All necessary commands for running SAIGE and PLINK are provided in the `.sh` files within this repository.

In addition to SAIGE and PLINK, the scripts in this repository require R and Python, along with the following libraries:

### R Libraries

Make sure the following R libraries are installed:

- `ggplot2`
- `dplyr`
- `tidyr`
- `gridExtra`
- `MASS`

You can install these libraries in R using the command:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "gridExtra", "MASS"))
```


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


## Data Preprocessing

Before running the analysis scripts, ensure that the raw data is processed correctly. Use the scripts in the `data_preprocess/` directory to clean and transform the data as needed.


## Single Variant Analysis

To conduct single variant analysis, navigate to the `Single_variant_analysis/linear_regression_classification/` directory and run the scripts provided.

# Run PLINK for single variant analysis
bash Single_variant_analysis/linear_regression_classification/plink_code.sh

## Combined Variants Analysis

Combined variant analysis includes the use of statistical tests like Burden, SKAT, and SKAT-O. Navigate to the `Combined_variants_analysis/` directory and run the respective scripts.


# Run SAIGE for combined variants analysis

## Contributing

Contributions are welcome! Please fork this repository and submit a pull request for any proposed changes or additions.

## License
Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
