# ============================================================
# 03_random_forest_importance.R
# Random forest metabolite importance ranking
# ============================================================

rank_dat <- combined_log2_df %>%
  dplyr::filter(Type %in% RF_IMPORTANCE_SET, Group %in% c(POS_CLASS, NEG_CLASS))

rank_dat$Group <- factor(rank_dat$Group, levels = c(POS_CLASS, NEG_CLASS))

X_rank <- rank_dat[, feature_cols, drop = FALSE]
y_rank <- rank_dat$Group

set.seed(123)
rf_rank_fit <- randomForest::randomForest(
  x = X_rank,
  y = y_rank,
  ntree = NTREE,
  importance = TRUE
)

imp <- randomForest::importance(rf_rank_fit)

importance_df <- data.frame(
  Rank = seq_len(nrow(imp)),
  Metabolite = rownames(imp),
  MeanDecreaseGini = imp[, "MeanDecreaseGini"],
  stringsAsFactors = FALSE,
  row.names = NULL
) %>%
  dplyr::left_join(
    metabolite_annotation %>% dplyr::select(Metabolite, Ion_Mode, feature, Metabolite_name),
    by = "Metabolite"
  ) %>%
  dplyr::arrange(dplyr::desc(MeanDecreaseGini)) %>%
  dplyr::mutate(Rank = dplyr::row_number()) %>%
  dplyr::select(Rank, Metabolite, Ion_Mode, feature, Metabolite_name, MeanDecreaseGini)

openxlsx::write.xlsx(
  importance_df,
  file.path(OUT_DIR, "02_RF_metabolite_importance_rankingall.xlsx"),
  asTable = TRUE
)

rf_rank <- importance_df$Metabolite

saveRDS(
  list(
    rf_rank_fit = rf_rank_fit,
    importance_df = importance_df,
    rf_rank = rf_rank
  ),
  file.path(OUT_DIR, "intermediate_03_rf_importance.rds")
)

log_msg("Random forest importance ranking completed.")
