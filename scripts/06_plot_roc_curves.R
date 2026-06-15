# ============================================================
# 06_plot_roc_curves.R
# Plot validation and test ROC curves
# ============================================================

plot_external_rocs <- function(score_df, dataset_type, out_pdf, out_png) {
  score_sub <- score_df %>% dplyr::filter(Type == dataset_type)
  
  if (nrow(score_sub) == 0) {
    warning("No data found for Type = ", dataset_type, "; skipping ROC plot.")
    return(NULL)
  }
  
  roc_plot_list <- list()
  auc_label_list <- list()
  
  for (k in sort(unique(score_sub$Feature_N))) {
    dat_k <- score_sub %>% dplyr::filter(Feature_N == k)
    
    rr <- safe_roc_auc_ci(
      response = dat_k$Group,
      predictor = dat_k$Score,
      positive_class = POS_CLASS,
      negative_class = NEG_CLASS
    )
    
    if (is.null(rr$roc)) next
    
    roc_plot_list[[as.character(k)]] <- data.frame(
      Dataset = dataset_type,
      Feature_N = as.character(k),
      FPR = 1 - rr$roc$specificities,
      TPR = rr$roc$sensitivities,
      AUC = rr$auc,
      CI_low = rr$ci_low,
      CI_high = rr$ci_high
    )
    
    auc_label_list[[as.character(k)]] <- data.frame(
      Dataset = dataset_type,
      Feature_N = k,
      AUC = rr$auc,
      CI_low = rr$ci_low,
      CI_high = rr$ci_high
    )
  }
  
  plot_df <- dplyr::bind_rows(roc_plot_list)
  label_df <- dplyr::bind_rows(auc_label_list)
  
  if (nrow(plot_df) == 0) {
    warning(dataset_type, " had no ROC curves to plot.")
    return(NULL)
  }
  
  label_text <- label_df %>%
    dplyr::arrange(Feature_N) %>%
    dplyr::mutate(
      Label = paste0(
        "Top ", Feature_N,
        ": AUC=", sprintf("%.3f", AUC),
        " (95% CI ", sprintf("%.3f", CI_low), "-", sprintf("%.3f", CI_high), ")"
      )
    ) %>%
    dplyr::pull(Label) %>%
    paste(collapse = "\n")
  
  p <- ggplot(plot_df, aes(x = FPR, y = TPR, color = Feature_N)) +
    geom_line(linewidth = 1.0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    annotate("text", x = 0.55, y = 0.25, label = label_text, hjust = 0, vjust = 1, size = 3.4, color = "black") +
    coord_equal() +
    labs(
      title = paste0(dataset_type, " ROC curves"),
      subtitle = "Train-fitted Random Forest models applied to external dataset",
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Top N"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")
  
  ggsave(out_pdf, p, width = 8.5, height = 6.5)
  ggsave(out_png, p, width = 8.5, height = 6.5, dpi = 600)
  
  return(p)
}

p_val_roc <- plot_external_rocs(
  score_df = val_test_scores_all,
  dataset_type = "Val",
  out_pdf = file.path(OUT_DIR, "07_Val_ROC_top5_10_15_25_100.pdf"),
  out_png = file.path(OUT_DIR, "07_Val_ROC_top5_10_15_25_100.png")
)

p_test_roc <- plot_external_rocs(
  score_df = val_test_scores_all,
  dataset_type = "Test",
  out_pdf = file.path(OUT_DIR, "08_Test_ROC_top5_10_15_25_100.pdf"),
  out_png = file.path(OUT_DIR, "08_Test_ROC_top5_10_15_25_100.png")
)
