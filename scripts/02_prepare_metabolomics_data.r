# ============================================================
# Prepare metabolomics matrix
# ============================================================

cat(
  "",
  file = file.path(
    OUT_DIR,
    "RF_metabolomics_run.log"
  )
)

sample_info_raw <-
  readxl::read_excel(
    SAMPLE_XLSX,
    sheet = SAMPLE_SHEET,
    col_types = "text"
  )

sample_info_raw <-
  as.data.frame(
    sample_info_raw,
    stringsAsFactors = FALSE
  )

no_col <- find_col(
  sample_info_raw,
  c("No.", "No", "no", "ID", "Person_ID")
)

pos_id_col <- find_col(
  sample_info_raw,
  c(
    "Sample_IDPOS",
    "SampleIDPOS",
    "POS_Sample_ID",
    "POS",
    "Sample_ID_POS"
  )
)

neg_id_col <- find_col(
  sample_info_raw,
  c(
    "Sample_IDNEG",
    "SampleIDNEG",
    "NEG_Sample_ID",
    "NEG",
    "Sample_ID_NEG"
  )
)

type_col <- find_col(
  sample_info_raw,
  c("Type", "type")
)

group_col <- find_col(
  sample_info_raw,
  c("Group", "group")
)

sample_info <- sample_info_raw %>%
  dplyr::transmute(
    No =
      normalize_id(
        .data[[no_col]]
      ),
    Sample_IDPOS =
      normalize_id(
        .data[[pos_id_col]]
      ),
    Sample_IDNEG =
      normalize_id(
        .data[[neg_id_col]]
      ),
    Type =
      normalize_id(
        .data[[type_col]]
      ),
    Group =
      normalize_id(
        .data[[group_col]]
      )
  ) %>%
  dplyr::filter(
    !is.na(No),
    No != ""
  ) %>%
  dplyr::distinct(
    No,
    .keep_all = TRUE
  )

if (
  anyDuplicated(
    sample_info$No
  ) > 0
) {
  stop(
    "Duplicated subject identifiers were found in sample_information.xlsx.",
    call. = FALSE
  )
}

if (
  !all(
    c(
      POS_CLASS,
      NEG_CLASS
    ) %in%
    unique(
      sample_info$Group
    )
  )
) {
  
  stop(
    "Both GC and NGC labels are required in the Group column.",
    call. = FALSE
  )
}

log_msg(
  "Sample information loaded: ",
  nrow(sample_info),
  " subjects"
)

log_msg(
  "Dataset distribution: ",
  paste(
    names(table(sample_info$Type)),
    as.integer(table(sample_info$Type)),
    collapse = "; "
  )
)

log_msg(
  "Group distribution: ",
  paste(
    names(table(sample_info$Group)),
    as.integer(table(sample_info$Group)),
    collapse = "; "
  )
)

pos_obj <- build_ion_matrix(
  POS_PEAK_CSV,
  POS_DIFF_CSV,
  "POS",
  sample_info,
  "Sample_IDPOS"
)

neg_obj <- build_ion_matrix(
  NEG_PEAK_CSV,
  NEG_DIFF_CSV,
  "NEG",
  sample_info,
  "Sample_IDNEG"
)

meta_cols <- c(
  "No",
  "Sample_IDPOS",
  "Sample_IDNEG",
  "Type",
  "Group"
)

pos_mat <- pos_obj$sample_matrix
neg_mat <- neg_obj$sample_matrix

metabolite_annotation <-
  dplyr::bind_rows(
    pos_obj$annotation,
    neg_obj$annotation
  ) %>%
  dplyr::distinct(
    Metabolite,
    .keep_all = TRUE
  )

combined_df <- sample_info %>%
  dplyr::select(
    dplyr::all_of(
      meta_cols
    )
  ) %>%
  dplyr::inner_join(
    pos_mat %>%
      dplyr::select(
        -Sample_IDPOS
      ),
    by = "No"
  ) %>%
  dplyr::inner_join(
    neg_mat %>%
      dplyr::select(
        -Sample_IDNEG
      ),
    by = "No"
  )

feature_cols <-
  setdiff(
    colnames(combined_df),
    meta_cols
  )

if (
  length(feature_cols) == 0
) {
  stop(
    "No metabolite features remained after merging POS and NEG data.",
    call. = FALSE
  )
}

for (fc in feature_cols) {
  
  combined_df[[fc]] <-
    suppressWarnings(
      as.numeric(
        combined_df[[fc]]
      )
    )
  
  x <- combined_df[[fc]]
  
  if (all(is.na(x))) {
    
    combined_df[[fc]] <- 0
    
  } else {
    
    min_non_na <- min(
      x[!is.na(x)],
      na.rm = TRUE
    )
    
    fill_value <- ifelse(
      is.finite(min_non_na),
      min_non_na / 2,
      0
    )
    
    combined_df[[fc]][
      is.na(
        combined_df[[fc]]
      )
    ] <- fill_value
  }
  
  combined_df[[fc]][
    combined_df[[fc]] < 0
  ] <- 0
}

combined_log2_df <- combined_df

combined_log2_df[
  feature_cols
] <- log2(
  combined_log2_df[
    feature_cols
  ] + 1
)

openxlsx::write.xlsx(
  combined_log2_df,
  file.path(
    OUT_DIR,
    "01_combined_diff_metabolite_log2_matrix.xlsx"
  ),
  asTable = TRUE
)

saveRDS(
  combined_log2_df,
  file.path(
    OUT_DIR,
    "combined_log2_matrix.rds"
  )
)

saveRDS(
  metabolite_annotation,
  file.path(
    OUT_DIR,
    "metabolite_annotation.rds"
  )
)

log_msg(
  "Combined samples = ",
  nrow(
    combined_log2_df
  )
)

log_msg(
  "Metabolite features = ",
  length(
    feature_cols
  )
)