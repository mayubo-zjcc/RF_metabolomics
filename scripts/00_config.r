# ============================================================
# Configuration
# ============================================================

required_packages <- c(
  "tidyverse",
  "readxl",
  "stringi",
  "openxlsx",
  "randomForest",
  "caret",
  "pROC",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages before running the pipeline: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(
  lapply(
    required_packages,
    function(pkg) {
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)
      )
    }
  )
)

set.seed(123)

# Run the pipeline from the repository root directory.
PROJECT_DIR <- normalizePath(".", mustWork = FALSE)

DATA_DIR <- file.path(PROJECT_DIR, "data")
OUT_DIR  <- file.path(PROJECT_DIR, "results")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Input files
POS_PEAK_CSV <- file.path(
  DATA_DIR,
  "rsd30posdata.csv"
)

NEG_PEAK_CSV <- file.path(
  DATA_DIR,
  "rsd30negdata.csv"
)

POS_DIFF_CSV <- file.path(
  DATA_DIR,
  "POS_differential_metabolites_q0.05.csv"
)

NEG_DIFF_CSV <- file.path(
  DATA_DIR,
  "NEG_differential_metabolites_q0.05.csv"
)

SAMPLE_XLSX <- file.path(
  DATA_DIR,
  "sample_information.xlsx"
)

SAMPLE_SHEET <- 1

# Classification labels
POS_CLASS <- "GC"
NEG_CLASS <- "NGC"

# Random forest parameters
FEATURE_SIZES <- c(5, 10, 15, 25, 100)
NTREE <- 1000
CV_FOLDS <- 5

# Dataset definitions
RF_IMPORTANCE_SET <- c("Train")
TRAIN_TYPE <- "Train"
VALID_TYPES <- c("Val", "Test")