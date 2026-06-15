# ============================================================
# 01_functions.R
# Helper functions
# ============================================================

normalize_id <- function(x) {
  x <- as.character(x)
  x %>%
    stringr::str_replace_all("\\u00A0", " ") %>%
    stringr::str_squish() %>%
    stringi::stri_trans_general("Fullwidth-Halfwidth")
}

make_safe_feature_name <- function(x, ion_prefix) {
  x <- normalize_id(x)
  x <- ifelse(is.na(x) | x == "", "UnknownFeature", x)
  x <- make.names(x, unique = FALSE)
  paste0(ion_prefix, "__", x)
}

find_col <- function(df, candidates, required = TRUE) {
  cn <- colnames(df)
  cn_norm <- tolower(normalize_id(cn))
  cand_norm <- tolower(normalize_id(candidates))
  idx <- match(cand_norm, cn_norm)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) > 0) {
    return(cn[idx[1]])
  }
  
  if (required) {
    stop("Missing required column: ", paste(candidates, collapse = " / "), call. = FALSE)
  }
  
  return(NA_character_)
}

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = RUN_LOG, append = TRUE)
}

read_peak_csv_auto_header <- function(csv_file, sample_ids, ion_label, max_header_scan = 20) {
  raw <- suppressWarnings(readr::read_csv(
    csv_file,
    col_names = FALSE,
    show_col_types = FALSE,
    guess_max = 100000,
    progress = FALSE
  ))
  
  if (nrow(raw) < 2) {
    stop(ion_label, ": the peak CSV has too few rows to read.", call. = FALSE)
  }
  
  sample_ids_norm <- normalize_id(sample_ids)
  n_scan <- min(max_header_scan, nrow(raw))
  
  scores <- purrr::map_dbl(seq_len(n_scan), function(i) {
    vals <- normalize_id(unlist(raw[i, ], use.names = FALSE))
    vals_lc <- tolower(vals)
    has_name <- any(vals_lc == "name" | vals_lc == "feature" | vals_lc == "metabolite")
    sample_match_n <- sum(vals %in% sample_ids_norm, na.rm = TRUE)
    sample_match_n + ifelse(has_name, 1000, 0)
  })
  
  header_row <- which.max(scores)
  header_vals <- normalize_id(unlist(raw[header_row, ], use.names = FALSE))
  header_vals[is.na(header_vals) | header_vals == ""] <- paste0("X", which(is.na(header_vals) | header_vals == ""))
  header_vals <- make.unique(header_vals)
  
  dat <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
  colnames(dat) <- header_vals
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  dat <- dat[rowSums(!is.na(dat) & dat != "") > 0, , drop = FALSE]
  
  log_msg(ion_label, ": detected peak table header row = ", header_row)
  log_msg(ion_label, ": number of matched sample columns in header = ", sum(header_vals %in% sample_ids_norm, na.rm = TRUE))
  
  return(dat)
}

read_diff_features <- function(diff_csv, ion_label) {
  diff_df <- suppressWarnings(readr::read_csv(
    diff_csv,
    show_col_types = FALSE,
    guess_max = 100000,
    progress = FALSE
  ))
  
  feature_col <- find_col(
    diff_df,
    c("feature", "Feature", "name", "Name", "metabolite", "Metabolite")
  )
  
  metabolite_name_col <- find_col(
    diff_df,
    c(
      "Metabolite_name", "metabolite_name",
      "Metabolite Name", "metabolite name",
      "MetaboliteName", "metaboliteName",
      "Name", "name",
      "Compound", "compound",
      "Compound_name", "compound_name",
      "Annotation", "annotation"
    ),
    required = FALSE
  )
  
  features <- normalize_id(diff_df[[feature_col]])
  keep <- !is.na(features) & features != ""
  
  diff_anno <- diff_df[keep, , drop = FALSE]
  diff_anno$feature <- normalize_id(diff_anno[[feature_col]])
  
  if (!is.na(metabolite_name_col)) {
    diff_anno$Metabolite_name <- normalize_id(diff_anno[[metabolite_name_col]])
  } else {
    diff_anno$Metabolite_name <- NA_character_
  }
  
  diff_anno <- diff_anno %>%
    dplyr::mutate(
      Ion_Mode = ion_label,
      Metabolite = make_safe_feature_name(feature, ion_label)
    ) %>%
    dplyr::select(Metabolite, Ion_Mode, feature, Metabolite_name, dplyr::everything()) %>%
    dplyr::distinct(Metabolite, .keep_all = TRUE)
  
  features <- unique(diff_anno$feature)
  log_msg(ion_label, ": number of differential metabolite features = ", length(features))
  
  return(list(features = features, annotation = diff_anno))
}

build_ion_matrix <- function(peak_csv, diff_csv, ion_label, sample_info, sample_col_name) {
  log_msg("\n==============================")
  log_msg("Reading ", ion_label, " data")
  log_msg("Peak file: ", peak_csv)
  log_msg("Differential metabolite file: ", diff_csv)
  log_msg("Sample column in sample information file: ", sample_col_name)
  
  sample_ids <- sample_info[[sample_col_name]]
  sample_ids_norm <- normalize_id(sample_ids)
  
  peak_df <- read_peak_csv_auto_header(peak_csv, sample_ids_norm, ion_label)
  diff_obj <- read_diff_features(diff_csv, ion_label)
  diff_features <- diff_obj$features
  diff_annotation <- diff_obj$annotation
  
  name_col <- find_col(peak_df, c("name", "Name", "feature", "Feature", "metabolite", "Metabolite"))
  
  peak_df$.feature_raw <- normalize_id(peak_df[[name_col]])
  peak_df$.feature_key <- tolower(peak_df$.feature_raw)
  diff_key <- tolower(normalize_id(diff_features))
  
  peak_sub <- peak_df %>% dplyr::filter(.feature_key %in% diff_key)
  
  if (nrow(peak_sub) == 0) {
    stop(ion_label, ": no metabolites were matched between the differential metabolite feature column and the peak table name column.", call. = FALSE)
  }
  
  log_msg(ion_label, ": matched differential metabolite rows = ", nrow(peak_sub))
  log_msg(ion_label, ": matched unique metabolite names = ", length(unique(peak_sub$.feature_raw)))
  
  peak_col_norm <- normalize_id(colnames(peak_sub))
  match_idx <- match(sample_ids_norm, peak_col_norm)
  
  matched_n <- sum(!is.na(match_idx))
  log_msg(ion_label, ": number of samples in sample information file = ", length(sample_ids_norm))
  log_msg(ion_label, ": number of matched sample columns in peak table = ", matched_n)
  
  if (matched_n == 0) {
    log_msg(ion_label, ": first 30 peak table column names: ", paste(head(colnames(peak_sub), 30), collapse = " | "))
    log_msg(ion_label, ": first 30 sample names: ", paste(head(sample_ids_norm, 30), collapse = " | "))
    stop(ion_label, ": peak table column names did not match ", sample_col_name, " in the sample information file.", call. = FALSE)
  }
  
  keep_sample_rows <- which(!is.na(match_idx))
  keep_peak_cols <- colnames(peak_sub)[match_idx[keep_sample_rows]]
  
  value_df <- peak_sub %>% dplyr::select(.feature_raw, dplyr::all_of(keep_peak_cols))
  
  for (cc in keep_peak_cols) {
    value_df[[cc]] <- suppressWarnings(as.numeric(value_df[[cc]]))
  }
  
  value_df <- value_df %>%
    dplyr::group_by(.feature_raw) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(keep_peak_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  feature_names <- make_safe_feature_name(value_df$.feature_raw, ion_label)
  mat <- as.matrix(value_df[, keep_peak_cols, drop = FALSE])
  rownames(mat) <- feature_names
  
  sample_mat <- as.data.frame(t(mat), stringsAsFactors = FALSE)
  sample_mat[] <- lapply(sample_mat, function(z) suppressWarnings(as.numeric(z)))
  
  sample_out <- sample_info[keep_sample_rows, c("No", sample_col_name), drop = FALSE]
  colnames(sample_out)[colnames(sample_out) == sample_col_name] <- paste0("Sample_ID", ion_label)
  sample_out <- dplyr::bind_cols(sample_out, sample_mat)
  
  return(list(sample_matrix = sample_out, annotation = diff_annotation))
}

safe_roc_auc_ci <- function(response, predictor, positive_class = POS_CLASS, negative_class = NEG_CLASS) {
  response <- factor(response, levels = c(positive_class, negative_class))
  
  if (length(unique(stats::na.omit(response))) < 2) {
    return(list(roc = NULL, auc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  
  roc_obj <- pROC::roc(
    response = response,
    predictor = predictor,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci_val <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95))
  
  return(list(roc = roc_obj, auc = auc_val, ci_low = ci_val[1], ci_high = ci_val[3]))
}

plot_train_cv_rocs <- function(roc_list, out_pdf, out_png) {
  plot_df <- purrr::imap_dfr(roc_list, function(x, nm) {
    if (is.null(x$roc)) return(NULL)
    data.frame(
      Feature_N = nm,
      FPR = 1 - x$roc$specificities,
      TPR = x$roc$sensitivities,
      AUC = x$auc,
      CI_low = x$ci_low,
      CI_high = x$ci_high
    )
  })
  
  label_df <- plot_df %>%
    dplyr::distinct(Feature_N, AUC, CI_low, CI_high) %>%
    dplyr::mutate(
      Label = paste0(
        "Top ", Feature_N,
        ": AUC=", sprintf("%.3f", AUC),
        " (95% CI ", sprintf("%.3f", CI_low), "-", sprintf("%.3f", CI_high), ")"
      )
    )
  
  label_text <- paste(label_df$Label, collapse = "\n")
  
  p <- ggplot(plot_df, aes(x = FPR, y = TPR, color = Feature_N)) +
    geom_line(linewidth = 1.0) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    annotate("text", x = 0.55, y = 0.25, label = label_text, hjust = 0, vjust = 1, size = 3.4, color = "black") +
    coord_equal() +
    labs(
      title = "Train 5-fold CV ROC curves",
      subtitle = "Random Forest using top ranked differential metabolites",
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

calc_binom_ci <- function(x, n, conf.level = 0.95) {
  if (is.na(x) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  ci <- stats::binom.test(x = x, n = n, conf.level = conf.level)$conf.int
  return(c(as.numeric(ci[1]), as.numeric(ci[2])))
}

calc_classification_metrics <- function(response, predictor, dataset_name, feature_n,
                                        positive_class = POS_CLASS,
                                        negative_class = NEG_CLASS,
                                        conf.level = 0.95) {
  response <- factor(response, levels = c(positive_class, negative_class))
  predictor <- as.numeric(predictor)
  
  if (length(unique(stats::na.omit(response))) < 2) {
    return(data.frame(
      Dataset = dataset_name,
      Feature_N = feature_n,
      Cutoff = NA_real_,
      AUC = NA_real_,
      AUC_CI_low = NA_real_,
      AUC_CI_high = NA_real_,
      Sensitivity = NA_real_,
      Sensitivity_CI_low = NA_real_,
      Sensitivity_CI_high = NA_real_,
      Specificity = NA_real_,
      Specificity_CI_low = NA_real_,
      Specificity_CI_high = NA_real_,
      Accuracy = NA_real_,
      Accuracy_CI_low = NA_real_,
      Accuracy_CI_high = NA_real_,
      TP = NA_integer_, TN = NA_integer_, FP = NA_integer_, FN = NA_integer_,
      N = length(response)
    ))
  }
  
  roc_obj <- pROC::roc(
    response = response,
    predictor = predictor,
    levels = c(negative_class, positive_class),
    direction = "<",
    quiet = TRUE
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- as.numeric(pROC::ci.auc(roc_obj, conf.level = conf.level))
  
  best_coords <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  cutoff <- as.numeric(best_coords[["threshold"]])
  pred_class <- ifelse(predictor >= cutoff, positive_class, negative_class)
  pred_class <- factor(pred_class, levels = c(positive_class, negative_class))
  
  TP <- sum(response == positive_class & pred_class == positive_class, na.rm = TRUE)
  TN <- sum(response == negative_class & pred_class == negative_class, na.rm = TRUE)
  FP <- sum(response == negative_class & pred_class == positive_class, na.rm = TRUE)
  FN <- sum(response == positive_class & pred_class == negative_class, na.rm = TRUE)
  N <- TP + TN + FP + FN
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  accuracy <- ifelse(N > 0, (TP + TN) / N, NA_real_)
  
  sensitivity_ci <- calc_binom_ci(TP, TP + FN, conf.level = conf.level)
  specificity_ci <- calc_binom_ci(TN, TN + FP, conf.level = conf.level)
  accuracy_ci <- calc_binom_ci(TP + TN, N, conf.level = conf.level)
  
  data.frame(
    Dataset = dataset_name,
    Feature_N = feature_n,
    Cutoff = cutoff,
    AUC = auc_val,
    AUC_CI_low = auc_ci[1],
    AUC_CI_high = auc_ci[3],
    Sensitivity = sensitivity,
    Sensitivity_CI_low = sensitivity_ci[1],
    Sensitivity_CI_high = sensitivity_ci[2],
    Specificity = specificity,
    Specificity_CI_low = specificity_ci[1],
    Specificity_CI_high = specificity_ci[2],
    Accuracy = accuracy,
    Accuracy_CI_low = accuracy_ci[1],
    Accuracy_CI_high = accuracy_ci[2],
    TP = TP, TN = TN, FP = FP, FN = FN, N = N
  )
}
