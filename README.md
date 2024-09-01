# Investigating the effects of coding *PCSK9* gene variants on blood glucose levels and susceptibility to type 2 diabetes

This repository contains scripts and data preprocessing workflows for analyzing the association between PCSK9 variants and Type 2 Diabetes (T2D) risk, including random glucose levels (RG) and lipid traits. The analysis is conducted using various statistical methods, including single-variant and combined-variant approaches.

## Table of Contents

- [Important Information Before Reproducing the Project](#important-information-before-reproducing-the-project)
- [Setup and Installation](#setup-and-installation)
- [Project Structure](#project-structure)
- [Data Preprocessing](#data-preprocessing)
- [Single Variant Analysis](#single-variant-analysis)
- [Combined Variants Analysis](#combined-variants-analysis)
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

- `deal_raw_data.sh`: Shell script to preprocess the raw data files.
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

- `hierarchical_clustering.R`: R script for hierarchical clustering of genetic variants based on Z-scores.
- `SAIGE_scripts.sh`: Shell script to run SAIGE for group-based analysis.
- `analysis_combined_pvalues.R`: R script to analyze combined p-values from multiple tests.

## Data Preprocessing

After downloading the UK Biobank data, you should have a `.vcf` zip file containing the genetic information for all individuals, as well as a separate file containing phenotypic data (including Type 2 diabetes status, BMI, LDLC, and random glucose levels).

To begin preprocessing, navigate to the `data_preprocess/` directory.

- Run the `deal_raw_data.sh` script using PLINK to process the raw genetic data.
- The other scripts in this directory can be easily run in R to handle additional preprocessing tasks.

Make sure you have the required software installed and configured correctly before running these scripts.

## Single Variant Analysis

To conduct a single variant analysis, navigate to the `Single_variant_analysis/` directory. This analysis consists of two parts:

1. **Linear Regression and Classification**
2. **Regression and Classification with L1 Penalty**

### 1. Linear Regression and Classification

To perform linear regression and classification:

- Navigate to the `linear_regression_classification/` directory.
- Run the `plink_code.sh` script using PLINK to perform the analysis.
- After running the PLINK script, use `significant_variants.R` to output the significant results.

### 2. Regression and Classification with L1 Penalty

To perform regression and classification with an L1 penalty:

- Navigate to the `regression_classification_penalty/` directory.
- Use Python to first encode all the data by running the `encode_data.py` script. Note that the data needed for `encode_data.py` should be prepared during the data preprocessing step.
- Run the other Python scripts in this directory to complete the analysis with the L1 penalty.


## Combined Variants Analysis

To perform a combined variants analysis, follow these steps:

1. **Navigate to the Analysis Directory**  
   Start by navigating to the `Combined_variants_analysis/` directory.

2. **Perform Hierarchical Clustering**  
   Run the `hierarchical_clustering.R` script to perform hierarchical clustering. This script will create a clustered file that will be used in the subsequent steps.

3. **Run SAIGE for Combined Variants Analysis**  
   Execute the `SAIGE_scripts.sh` script to run SAIGE for the combined variants analysis. Run each line of the script sequentially using the SAIGE software.

4. **Analyze Combined P-values**  
   After running the SAIGE analysis, use the `analysis_combined_pvalues.R` script to analyze the combined p-values and interpret the results.

### Notes

- Ensure that the SAIGE software is properly installed and configured before running the analysis.
- Make sure all necessary input files are prepared and available in the appropriate directories before starting the analysis.


## License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

