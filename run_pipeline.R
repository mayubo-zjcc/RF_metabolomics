# ============================================================
# Run the complete random forest metabolomics pipeline
# ============================================================

source("scripts/00_config.R")
source("scripts/01_functions.R")
source("scripts/02_prepare_metabolomics_matrix.R")
source("scripts/03_random_forest_importance.R")
source("scripts/04_train_cross_validation.R")
source("scripts/05_score_validation_test.R")
source("scripts/06_plot_roc_curves.R")
source("scripts/07_export_performance_metrics.R")

message("Random forest metabolomics pipeline finished successfully.")
