train_dat <- combined_log2_df %>%
  dplyr::filter(
    Type == TRAIN_TYPE,
    Group %in%
      c(
        POS_CLASS,
        NEG_CLASS
      )
  )

train_dat$Group <- factor(
  train_dat$Group,
  levels = c(
    POS_CLASS,
    NEG_CLASS
  )
)

feature_sizes <-
  FEATURE_SIZES[
    FEATURE_SIZES <=
      length(rf_rank)
  ]

ctrl_cv <- caret::trainControl(
  method = "cv",
  number = CV_FOLDS,
  classProbs = TRUE,
  summaryFunction =
    caret::twoClassSummary,
  savePredictions = "final",
  allowParallel = TRUE
)

cv_scores_list <- list()
cv_auc_list <- list()
roc_list <- list()
feature_list <- list()

for (k in feature_sizes) {
  
  log_msg(
    "Train 5-fold CV: top ",
    k,
    " metabolites"
  )
  
  top_features <-
    rf_rank[1:k]
  
  feature_list[
    [paste0("Top", k)]
  ] <- data.frame(
    Rank =
      seq_along(
        top_features
      ),
    Metabolite =
      top_features
  ) %>%
    dplyr::left_join(
      metabolite_annotation %>%
        dplyr::select(
          Metabolite,
          Ion_Mode,
          feature,
          Metabolite_name
        ),
      by = "Metabolite"
    )
  
  X_train <-
    train_dat[
      ,
      top_features,
      drop = FALSE
    ]
  
  y_train <-
    train_dat$Group
  
  tune_grid <-
    expand.grid(
      mtry =
        max(
          1,
          floor(
            sqrt(k)
          )
        )
    )
  
  set.seed(
    123 + k
  )
  
  rf_cv_fit <-
    caret::train(
      x = X_train,
      y = y_train,
      method = "rf",
      metric = "ROC",
      trControl = ctrl_cv,
      tuneGrid = tune_grid,
      ntree = NTREE,
      importance = TRUE
    )
  
  pred_df <-
    rf_cv_fit$pred %>%
    dplyr::filter(
      mtry ==
        tune_grid$mtry
    ) %>%
    dplyr::mutate(
      Feature_N = k,
      Score =
        .data[[POS_CLASS]]
    )
  
  cv_score_mean <-
    pred_df %>%
    dplyr::group_by(
      rowIndex
    ) %>%
    dplyr::summarise(
      No =
        train_dat$No[
          first(rowIndex)
        ],
      Sample_IDPOS =
        train_dat$Sample_IDPOS[
          first(rowIndex)
        ],
      Sample_IDNEG =
        train_dat$Sample_IDNEG[
          first(rowIndex)
        ],
      Group =
        first(obs),
      Feature_N =
        first(Feature_N),
      CV_Score =
        mean(
          Score,
          na.rm = TRUE
        ),
      N_Pred =
        dplyr::n(),
      .groups = "drop"
    )
  
  roc_res <-
    safe_roc_auc_ci(
      response =
        cv_score_mean$Group,
      predictor =
        cv_score_mean$CV_Score
    )
  
  cv_auc_list[
    [as.character(k)]
  ] <- data.frame(
    Dataset = "Train_5fold_CV",
    Feature_N = k,
    AUC = roc_res$auc,
    CI_low = roc_res$ci_low,
    CI_high = roc_res$ci_high,
    mtry = tune_grid$mtry,
    ntree = NTREE
  )
  
  cv_scores_list[
    [as.character(k)]
  ] <- cv_score_mean
  
  roc_list[
    [as.character(k)]
  ] <- roc_res
}

cv_scores_all <-
  dplyr::bind_rows(
    cv_scores_list
  )

cv_auc_all <-
  dplyr::bind_rows(
    cv_auc_list
  )

openxlsx::write.xlsx(
  list(
    CV_scores =
      cv_scores_all,
    CV_AUC =
      cv_auc_all
  ),
  file.path(
    OUT_DIR,
    "04_Train_CV_scores_and_AUC.xlsx"
  ),
  asTable = TRUE
)

openxlsx::write.xlsx(
  feature_list,
  file.path(
    OUT_DIR,
    "06_topK_feature_lists.xlsx"
  ),
  asTable = TRUE
)