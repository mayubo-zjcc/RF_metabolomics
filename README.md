# Random forest metabolomics model pipeline

This repository contains R scripts used to preprocess POS/NEG metabolomics data, select differential metabolites, train random forest classifiers, evaluate 5-fold cross-validation performance, and score validation/test cohorts.

## Repository structure

# Random forest metabolomics model pipeline

This repository contains R scripts used to preprocess POS/NEG metabolomics data, select differential metabolites, train random forest classifiers, evaluate 5-fold cross-validation performance, and validate the models in two independent validation cohorts.

## Repository structure

```text
data/
  sample_information.xlsx
  rsd30posdata.csv
  rsd30negdata.csv
  POS_differential_metabolites_q0.05.csv
  NEG_differential_metabolites_q0.05.csv

scripts/
  00_config.R
  01_functions.R
  02_prepare_metabolomics_data.R
  03_random_forest_feature_ranking.R
  04_random_forest_cross_validation.R
  05_validation_analysis.R
  06_performance_metrics.R

results/
  Analysis outputs

## Required input files

Input files are stored in the data/ directory.
The sample information table should include:
- No: de-identified participant identifier
- Sample_IDPOS: sample identifier in the POS peak table
- Sample_IDNEG: sample identifier in the NEG peak table
- Type: dataset label (Train, Val, or Test)
- Group: class label (GC or NGC)
Val and Test represent two independent validation cohorts.

## How to run

Run the scripts in numerical order from the repository root directory:
source("scripts/00_config.R")
source("scripts/01_functions.R")
source("scripts/02_prepare_metabolomics_data.R")
source("scripts/03_random_forest_feature_ranking.R")
source("scripts/04_random_forest_cross_validation.R")
source("scripts/05_validation_analysis.R")
source("scripts/06_performance_metrics.R")

## Main outputs

Outputs are saved in the results/ directory.
- 01_combined_diff_metabolite_log2_matrix.xlsx
- 02_RF_metabolite_importance_rankingall.xlsx
- 03_Train_CV_ROC_top5_10_15_25_100.pdf/png
- 04_Train_CV_scores_and_AUC.xlsx
- 05_Val_Test_scores_top5_10_15_25_100.xlsx
- 06_topK_feature_lists.xlsx
- 07_Val_ROC_top5_10_15_25_100.pdf/png
- 08_Test_ROC_top5_10_15_25_100.pdf/png
- 10_Train_Val_Test_Performance_Metrics.xlsx

## R packages

The pipeline requires the following R packages:

- tidyverse
- readxl
- stringi
- openxlsx
- randomForest
- caret
- pROC
- ggplot2

Package installation is not performed automatically. Please install missing packages before running the pipeline.
