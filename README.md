# Random forest metabolomics model pipeline

This repository contains R scripts used to preprocess POS/NEG metabolomics data, select differential metabolites, train random forest classifiers, evaluate 5-fold cross-validation performance, and score validation/test cohorts.

## Repository structure

```text
scripts/
  00_config.R                         # File paths and analysis parameters
  01_functions.R                      # Helper functions
  02_prepare_metabolomics_matrix.R     # Read sample information and build POS/NEG log2 matrix
  03_random_forest_importance.R        # Train RF model for metabolite importance ranking
  04_train_cross_validation.R          # 5-fold cross-validation in the training set
  05_score_validation_test.R           # Train final RF models and score validation/test sets
  06_plot_roc_curves.R                 # Plot train/validation/test ROC curves
  07_export_performance_metrics.R      # Export AUC, sensitivity, specificity, and accuracy
run_pipeline.R                         # Run the complete workflow
```

## Required input files

Set the input file names in `scripts/00_config.R`.

Expected files:

1. POS peak intensity table
2. NEG peak intensity table
3. POS differential metabolite list
4. NEG differential metabolite list
5. Sample information table containing de-identified sample IDs, group labels, and dataset type labels

The sample information table should include columns equivalent to:

- `No`: participant/person identifier
- `Sample_IDPOS`: sample name in the POS peak table
- `Sample_IDNEG`: sample name in the NEG peak table
- `Type`: dataset split, such as `Train`, `Val`, or `Test`
- `Group`: class label, such as `GC` or `NGC`

## How to run

Edit `scripts/00_config.R`, or set the base directory before running:

```r
Sys.setenv(BASE_DIR = "/path/to/project")
source("run_pipeline.R")
```

## Main outputs

Outputs are saved to `OUT_DIR`, defined in `scripts/00_config.R`.

- `01_combined_diff_metabolite_log2_matrix.xlsx`
- `02_RF_metabolite_importance_rankingall.xlsx`
- `03_Train_CV_ROC_top5_10_15_25_100.pdf/png`
- `04_Train_CV_scores_and_AUC.xlsx`
- `05_Val_Test_scores_top5_10_15_25_100.xlsx`
- `06_topK_feature_lists.xlsx`
- `07_Val_ROC_top5_10_15_25_100.pdf/png`
- `08_Test_ROC_top5_10_15_25_100.pdf/png`
- `10_Train_Val_Test_Performance_Metrics.xlsx`

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
