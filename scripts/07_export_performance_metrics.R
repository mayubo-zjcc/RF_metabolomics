# ============================================================
# 07_export_performance_metrics.R
# Export sensitivity, specificity, accuracy, and AUC for Train/Val/Test
# ============================================================

train_metrics <- cv_scores_all %>%
  dplyr::group_split(Feature_N) %>%
  purrr::map_dfr(function(df) {
    calc_classification_metrics(
      response = df$Group,
      predictor = df$CV_Score,
      dataset_name = "Train_5fold_CV",
      feature_n = unique(df$Feature_N),
      positive_class = POS_CLASS,
      negative_class = NEG_CLASS
    )
  })

val_test_metrics <- val_test_scores_all %>%
  dplyr::group_split(Type, Feature_N) %>%
  purrr::map_dfr(function(df) {
    calc_classification_metrics(
      response = df$Group,
      predictor = df$Score,
      dataset_name = unique(df$Type),
      feature_n = unique(df$Feature_N),
      positive_class = POS_CLASS,
      negative_class = NEG_CLASS
    )
  })

all_performance_metrics <- dplyr::bind_rows(train_metrics, val_test_metrics) %>%
  dplyr::arrange(factor(Dataset, levels = c("Train_5fold_CV", "Val", "Test")), Feature_N)

all_performance_metrics_wide <- all_performance_metrics %>%
  dplyr::select(
    Dataset, Feature_N, Cutoff,
    AUC, AUC_CI_low, AUC_CI_high,
    Sensitivity, Sensitivity_CI_low, Sensitivity_CI_high,
    Specificity, Specificity_CI_low, Specificity_CI_high,
    Accuracy, Accuracy_CI_low, Accuracy_CI_high,
    TP, TN, FP, FN, N
  ) %>%
  tidyr::pivot_wider(
    names_from = Dataset,
    values_from = c(
      Cutoff,
      AUC, AUC_CI_low, AUC_CI_high,
      Sensitivity, Sensitivity_CI_low, Sensitivity_CI_high,
      Specificity, Specificity_CI_low, Specificity_CI_high,
      Accuracy, Accuracy_CI_low, Accuracy_CI_high,
      TP, TN, FP, FN, N
    )
  )

openxlsx::write.xlsx(
  list(
    Metrics_long = all_performance_metrics,
    Metrics_wide = all_performance_metrics_wide
  ),
  file.path(OUT_DIR, "10_Train_Val_Test_Performance_Metrics.xlsx"),
  asTable = TRUE
)

saveRDS(
  list(
    all_performance_metrics = all_performance_metrics,
    all_performance_metrics_wide = all_performance_metrics_wide
  ),
  file.path(OUT_DIR, "intermediate_07_performance_metrics.rds")
)

message("Done: Train/Val/Test sensitivity, specificity, accuracy, and AUC have been exported.")
message(file.path(OUT_DIR, "10_Train_Val_Test_Performance_Metrics.xlsx"))
