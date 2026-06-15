# ============================================================
# 00_config.R
# File paths and analysis parameters
# ============================================================

required_packages <- c(
  "tidyverse", "readxl", "stringi", "openxlsx",
  "randomForest", "caret", "pROC", "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages before running this pipeline: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

set.seed(123)

# Set BASE_DIR by either editing this line or using Sys.setenv(BASE_DIR = "/path/to/project").
BASE_DIR <- Sys.getenv("BASE_DIR", unset = ".")

POS_PEAK_CSV <- file.path(BASE_DIR, "rsd30posdata.csv")
NEG_PEAK_CSV <- file.path(BASE_DIR, "rsd30negdata.csv")

# Rename these folders/files to English when preparing a public GitHub repository.
POS_DIFF_CSV <- file.path(BASE_DIR, "differential_metabolites", "POS_differential_metabolites_q0.05.csv")
NEG_DIFF_CSV <- file.path(BASE_DIR, "differential_metabolites", "NEG_differential_metabolites_q0.05.csv")

SAMPLE_XLSX <- file.path(BASE_DIR, "Pan_Inf.xlsx")
SAMPLE_SHEET <- 1

OUT_DIR <- file.path(BASE_DIR, "outputs_RF_metabolomics")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

POS_CLASS <- "GC"
NEG_CLASS <- "NGC"

FEATURE_SIZES <- c(5, 10, 15, 25, 100)
NTREE <- 1000
CV_FOLDS <- 5

RF_IMPORTANCE_SET <- c("Train")
TRAIN_TYPE <- "Train"
VALID_TYPES <- c("Val", "Test")

RUN_LOG <- file.path(OUT_DIR, "RF_metabolomics_run.log")
cat("", file = RUN_LOG)
