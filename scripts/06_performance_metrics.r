# ============================================================
# Performance metrics
# ============================================================

calc_binom_ci <- function(
    x,
    n,
    conf.level = 0.95) {
  
  if (
    is.na(x) ||
    is.na(n) ||
    n <= 0
  ) {
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }
  
  ci <-
    stats::binom.test(
      x = x,
      n = n,
      conf.level = conf.level
    )$conf.int
  
  c(
    as.numeric(ci[1]),
    as.numeric(ci[2])
  )
}


calc_classification_metrics <- function(
    response,
    predictor,
    dataset_name,
    feature_n,
    positive_class = POS_CLASS,
    negative_class = NEG_CLASS,
    conf.level = 0.95) {
  
  response <- factor(
    response,
    levels = c(
      positive_class,
      negative_class
    )
  )
  
  predictor <-
    as.numeric(
      predictor
    )
  
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
  
  auc_ci <-
    as.numeric(
      pROC::ci.auc(
        roc_obj,
        conf.level = conf.level
      )
    )
  
  best_coords <-
    pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = c(
        "threshold",
        "sensitivity",
        "specificity"
      ),
      transpose = FALSE
    )
  
  cutoff <-
    as.numeric(
      best_coords[
        ["threshold"]
      ]
    )
  
  pred_class <- ifelse(
    predictor >= cutoff,
    positive_class,
    negative_class
  )
  
  pred_class <- factor(
    pred_class,
    levels = c(
      positive_class,
      negative_class
    )
  )
  
  TP <- sum(
    response == positive_class &
      pred_class == positive_class,
    na.rm = TRUE
  )
  
  TN <- sum(
    response == negative_class &
      pred_class == negative_class,
    na.rm = TRUE
  )
  
  FP <- sum(
    response == negative_class &
      pred_class == positive_class,
    na.rm = TRUE
  )
  
  FN <- sum(
    response == positive_class &
      pred_class == negative_class,
    na.rm = TRUE
  )
  
  N <- TP + TN + FP + FN
  
  sensitivity <-
    TP / (TP + FN)
  
  specificity <-
    TN / (TN + FP)
  
  accuracy <-
    (TP + TN) / N
  
  sensitivity_ci <-
    calc_binom_ci(
      TP,
      TP + FN,
      conf.level
    )
  
  specificity_ci <-
    calc_binom_ci(
      TN,
      TN + FP,
      conf.level
    )
  
  accuracy_ci <-
    calc_binom_ci(
      TP + TN,
      N,
      conf.level
    )
  
  data.frame(
    Dataset = dataset_name,
    Feature_N = feature_n,
    Cutoff = cutoff,
    
    AUC = auc_val,
    AUC_CI_low = auc_ci[1],
    AUC_CI_high = auc_ci[3],
    
    Sensitivity = sensitivity,
    Sensitivity_CI_low =
      sensitivity_ci[1],
    Sensitivity_CI_high =
      sensitivity_ci[2],
    
    Specificity = specificity,
    Specificity_CI_low =
      specificity_ci[1],
    Specificity_CI_high =
      specificity_ci[2],
    
    Accuracy = accuracy,
    Accuracy_CI_low =
      accuracy_ci[1],
    Accuracy_CI_high =
      accuracy_ci[2],
    
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    N = N
  )
}

train_metrics <-
  cv_scores_all %>%
  dplyr::group_split(
    Feature_N
  ) %>%
  purrr::map_dfr(
    function(df) {
      
      calc_classification_metrics(
        response =
          df$Group,
        predictor =
          df$CV_Score,
        dataset_name =
          "Train_5fold_CV",
        feature_n =
          unique(
            df$Feature_N
          )
      )
    }
  )

val_test_metrics <-
  val_test_scores_all %>%
  dplyr::group_split(
    Type,
    Feature_N
  ) %>%
  purrr::map_dfr(
    function(df) {
      
      calc_classification_metrics(
        response =
          df$Group,
        predictor =
          df$Score,
        dataset_name =
          unique(
            df$Type
          ),
        feature_n =
          unique(
            df$Feature_N
          )
      )
    }
  )

all_performance_metrics <-
  dplyr::bind_rows(
    train_metrics,
    val_test_metrics
  ) %>%
  dplyr::arrange(
    factor(
      Dataset,
      levels = c(
        "Train_5fold_CV",
        "Val",
        "Test"
      )
    ),
    Feature_N
  )

openxlsx::write.xlsx(
  list(
    Metrics_long =
      all_performance_metrics
  ),
  file.path(
    OUT_DIR,
    "10_Train_Val_Test_Performance_Metrics.xlsx"
  ),
  asTable = TRUE
)

message(
  "Performance metrics exported successfully."
)