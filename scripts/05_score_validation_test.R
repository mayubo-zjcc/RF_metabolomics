# ============================================================
# 05_score_validation_test.R
# Train final RF models and score validation/test sets
# ============================================================

valid_dat <- combined_log2_df %>%
  dplyr::filter(Type %in% VALID_TYPES, Group %in% c(POS_CLASS, NEG_CLASS))

valid_dat$Group <- factor(valid_dat$Group, levels = c(POS_CLASS, NEG_CLASS))

val_test_scores_list <- list()
val_test_auc_list <- list()

for (k in feature_sizes) {
  log_msg("\nTrain final RF and score Val/Test: top ", k, " metabolites")
  
  top_features <- rf_rank[1:k]
  
  X_train <- train_dat[, top_features, drop = FALSE]
  y_train <- train_dat$Group
  X_valid <- valid_dat[, top_features, drop = FALSE]
  
  set.seed(2024 + k)
  rf_fit <- randomForest::randomForest(
    x = X_train,
    y = y_train,
    ntree = NTREE,
    mtry = max(1, floor(sqrt(k))),
    importance = TRUE
  )
  
  valid_prob <- predict(rf_fit, X_valid, type = "prob")[, POS_CLASS]
  
  score_df <- valid_dat %>%
    dplyr::select(No, Sample_IDPOS, Sample_IDNEG, Type, Group) %>%
    dplyr::mutate(Feature_N = k, Score = as.numeric(valid_prob)) %>%
    dplyr::arrange(Type, dplyr::desc(Score))
  
  val_test_scores_list[[as.character(k)]] <- score_df
  
  auc_by_type <- score_df %>%
    dplyr::group_by(Type) %>%
    dplyr::group_modify(function(.x, .y) {
      rr <- safe_roc_auc_ci(
        response = .x$Group,
        predictor = .x$Score,
        positive_class = POS_CLASS,
        negative_class = NEG_CLASS
      )
      data.frame(
        Feature_N = k,
        AUC = rr$auc,
        CI_low = rr$ci_low,
        CI_high = rr$ci_high,
        ntree = NTREE,
        mtry = max(1, floor(sqrt(k)))
      )
    }) %>%
    dplyr::ungroup()
  
  val_test_auc_list[[as.character(k)]] <- auc_by_type
}

val_test_scores_all <- dplyr::bind_rows(val_test_scores_list)
val_test_auc_all <- dplyr::bind_rows(val_test_auc_list)

val_test_scores_wide <- val_test_scores_all %>%
  dplyr::select(No, Sample_IDPOS, Sample_IDNEG, Type, Group, Feature_N, Score) %>%
  tidyr::pivot_wider(names_from = Feature_N, values_from = Score, names_prefix = "Score_Top") %>%
  dplyr::arrange(Type, No)

openxlsx::write.xlsx(
  list(
    Val_Test_scores_long = val_test_scores_all,
    Val_Test_scores_wide = val_test_scores_wide,
    Val_Test_AUC = val_test_auc_all
  ),
  file.path(OUT_DIR, "05_Val_Test_scores_top5_10_15_25_100.xlsx"),
  asTable = TRUE
)

saveRDS(
  list(
    valid_dat = valid_dat,
    val_test_scores_all = val_test_scores_all,
    val_test_auc_all = val_test_auc_all,
    val_test_scores_wide = val_test_scores_wide
  ),
  file.path(OUT_DIR, "intermediate_05_val_test_scores.rds")
)

log_msg("Validation/test scoring completed.")
