# ============================================================
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
  x <- ifelse(
    is.na(x) | x == "",
    "UnknownFeature",
    x
  )
  
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
    stop(
      "Missing required column: ",
      paste(candidates, collapse = " / "),
      call. = FALSE
    )
  }
  
  return(NA_character_)
}

log_msg <- function(...) {
  
  msg <- paste0(...)
  
  cat(msg, "\n")
  
  cat(
    msg,
    "\n",
    file = file.path(
      OUT_DIR,
      "RF_metabolomics_run.log"
    ),
    append = TRUE
  )
}


read_peak_csv_auto_header <- function(
    csv_file,
    sample_ids,
    ion_label,
    max_header_scan = 20) {
  
  raw <- suppressWarnings(
    readr::read_csv(
      csv_file,
      col_names = FALSE,
      show_col_types = FALSE,
      guess_max = 100000,
      progress = FALSE
    )
  )
  
  if (nrow(raw) < 2) {
    stop(
      ion_label,
      ": peak table contains too few rows.",
      call. = FALSE
    )
  }
  
  sample_ids_norm <- normalize_id(sample_ids)
  
  n_scan <- min(
    max_header_scan,
    nrow(raw)
  )
  
  scores <- purrr::map_dbl(
    seq_len(n_scan),
    function(i) {
      
      vals <- normalize_id(
        unlist(raw[i, ], use.names = FALSE)
      )
      
      vals_lc <- tolower(vals)
      
      has_name <- any(
        vals_lc == "name" |
          vals_lc == "feature" |
          vals_lc == "metabolite"
      )
      
      sample_match_n <- sum(
        vals %in% sample_ids_norm,
        na.rm = TRUE
      )
      
      sample_match_n +
        ifelse(has_name, 1000, 0)
    }
  )
  
  header_row <- which.max(scores)
  
  header_vals <- normalize_id(
    unlist(
      raw[header_row, ],
      use.names = FALSE
    )
  )
  
  empty_idx <- which(
    is.na(header_vals) |
      header_vals == ""
  )
  
  header_vals[empty_idx] <- paste0(
    "X",
    empty_idx
  )
  
  header_vals <- make.unique(header_vals)
  
  dat <- raw[
    (header_row + 1):nrow(raw),
    ,
    drop = FALSE
  ]
  
  colnames(dat) <- header_vals
  
  dat <- as.data.frame(
    dat,
    stringsAsFactors = FALSE
  )
  
  dat <- dat[
    rowSums(
      !is.na(dat) &
        dat != ""
    ) > 0,
    ,
    drop = FALSE
  ]
  
  log_msg(
    ion_label,
    ": detected peak-table header row = ",
    header_row
  )
  
  log_msg(
    ion_label,
    ": matched sample columns = ",
    sum(
      header_vals %in% sample_ids_norm,
      na.rm = TRUE
    )
  )
  
  return(dat)
}


read_diff_features <- function(
    diff_csv,
    ion_label) {
  
  diff_df <- suppressWarnings(
    readr::read_csv(
      diff_csv,
      show_col_types = FALSE,
      guess_max = 100000,
      progress = FALSE
    )
  )
  
  feature_col <- find_col(
    diff_df,
    c(
      "feature",
      "Feature",
      "name",
      "Name",
      "metabolite",
      "Metabolite"
    )
  )
  
  metabolite_name_col <- find_col(
    diff_df,
    c(
      "Metabolite_name",
      "metabolite_name",
      "Metabolite Name",
      "metabolite name",
      "MetaboliteName",
      "metaboliteName",
      "Name",
      "name",
      "Compound",
      "compound",
      "Compound_name",
      "compound_name",
      "Annotation",
      "annotation"
    ),
    required = FALSE
  )
  
  features <- normalize_id(
    diff_df[[feature_col]]
  )
  
  keep <- !is.na(features) &
    features != ""
  
  diff_anno <- diff_df[
    keep,
    ,
    drop = FALSE
  ]
  
  diff_anno$feature <- normalize_id(
    diff_anno[[feature_col]]
  )
  
  if (!is.na(metabolite_name_col)) {
    
    diff_anno$Metabolite_name <-
      normalize_id(
        diff_anno[[metabolite_name_col]]
      )
    
  } else {
    
    diff_anno$Metabolite_name <-
      NA_character_
  }
  
  diff_anno <- diff_anno %>%
    dplyr::mutate(
      Ion_Mode = ion_label,
      Metabolite =
        make_safe_feature_name(
          feature,
          ion_label
        )
    ) %>%
    dplyr::select(
      Metabolite,
      Ion_Mode,
      feature,
      Metabolite_name,
      dplyr::everything()
    ) %>%
    dplyr::distinct(
      Metabolite,
      .keep_all = TRUE
    )
  
  features <- unique(
    diff_anno$feature
  )
  
  log_msg(
    ion_label,
    ": number of differential metabolite features = ",
    length(features)
  )
  
  return(
    list(
      features = features,
      annotation = diff_anno
    )
  )
}


build_ion_matrix <- function(
    peak_csv,
    diff_csv,
    ion_label,
    sample_info,
    sample_col_name) {
  
  log_msg(
    "\n=============================="
  )
  
  log_msg(
    "Processing ",
    ion_label,
    " data"
  )
  
  sample_ids <- sample_info[
    [sample_col_name]
  ]
  
  sample_ids_norm <- normalize_id(
    sample_ids
  )
  
  peak_df <- read_peak_csv_auto_header(
    peak_csv,
    sample_ids_norm,
    ion_label
  )
  
  diff_obj <- read_diff_features(
    diff_csv,
    ion_label
  )
  
  diff_features <- diff_obj$features
  diff_annotation <- diff_obj$annotation
  
  name_col <- find_col(
    peak_df,
    c(
      "name",
      "Name",
      "feature",
      "Feature",
      "metabolite",
      "Metabolite"
    )
  )
  
  peak_df$.feature_raw <-
    normalize_id(
      peak_df[[name_col]]
    )
  
  peak_df$.feature_key <-
    tolower(
      peak_df$.feature_raw
    )
  
  diff_key <- tolower(
    normalize_id(
      diff_features
    )
  )
  
  peak_sub <- peak_df %>%
    dplyr::filter(
      .feature_key %in% diff_key
    )
  
  if (nrow(peak_sub) == 0) {
    
    stop(
      ion_label,
      ": no differential metabolite features were matched to the peak table.",
      call. = FALSE
    )
  }
  
  peak_col_norm <- normalize_id(
    colnames(peak_sub)
  )
  
  match_idx <- match(
    sample_ids_norm,
    peak_col_norm
  )
  
  matched_n <- sum(
    !is.na(match_idx)
  )
  
  log_msg(
    ion_label,
    ": matched samples = ",
    matched_n
  )
  
  if (matched_n == 0) {
    
    stop(
      ion_label,
      ": sample identifiers could not be matched between sample_information.xlsx and the peak table.",
      call. = FALSE
    )
  }
  
  keep_sample_rows <- which(
    !is.na(match_idx)
  )
  
  keep_peak_cols <- colnames(
    peak_sub
  )[match_idx[keep_sample_rows]]
  
  value_df <- peak_sub %>%
    dplyr::select(
      .feature_raw,
      dplyr::all_of(
        keep_peak_cols
      )
    )
  
  for (cc in keep_peak_cols) {
    
    value_df[[cc]] <-
      suppressWarnings(
        as.numeric(
          value_df[[cc]]
        )
      )
  }
  
  value_df <- value_df %>%
    dplyr::group_by(
      .feature_raw
    ) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(
          keep_peak_cols
        ),
        ~ mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  feature_names <-
    make_safe_feature_name(
      value_df$.feature_raw,
      ion_label
    )
  
  mat <- as.matrix(
    value_df[
      ,
      keep_peak_cols,
      drop = FALSE
    ]
  )
  
  rownames(mat) <- feature_names
  
  sample_mat <- as.data.frame(
    t(mat),
    stringsAsFactors = FALSE
  )
  
  sample_mat[] <- lapply(
    sample_mat,
    function(z) {
      suppressWarnings(
        as.numeric(z)
      )
    }
  )
  
  sample_out <- sample_info[
    keep_sample_rows,
    c(
      "No",
      sample_col_name
    ),
    drop = FALSE
  ]
  
  colnames(sample_out)[
    colnames(sample_out) ==
      sample_col_name
  ] <- paste0(
    "Sample_ID",
    ion_label
  )
  
  sample_out <- dplyr::bind_cols(
    sample_out,
    sample_mat
  )
  
  return(
    list(
      sample_matrix = sample_out,
      annotation = diff_annotation
    )
  )
}


safe_roc_auc_ci <- function(
    response,
    predictor,
    positive_class = POS_CLASS,
    negative_class = NEG_CLASS) {
  
  response <- factor(
    response,
    levels = c(
      positive_class,
      negative_class
    )
  )
  
  if (
    length(
      unique(
        stats::na.omit(
          response
        )
      )
    ) < 2
  ) {
    
    return(
      list(
        roc = NULL,
        auc = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    )
  }
  
  roc_obj <- pROC::roc(
    response = response,
    predictor = predictor,
    levels = c(
      negative_class,
      positive_class
    ),
    direction = "<",
    quiet = TRUE
  )
  
  auc_val <-
    as.numeric(
      pROC::auc(
        roc_obj
      )
    )
  
  ci_val <-
    as.numeric(
      pROC::ci.auc(
        roc_obj,
        conf.level = 0.95
      )
    )
  
  list(
    roc = roc_obj,
    auc = auc_val,
    ci_low = ci_val[1],
    ci_high = ci_val[3]
  )
}