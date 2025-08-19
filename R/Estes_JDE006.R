#' Fit a Univariate Linear Model of CopyNumber Against a Single Feature
#' 
#' Customized for Estes SIV/HIV rebound study JDE006
#'
#' This function automates the process of fitting a simple linear regression
#' model where the response variable is `CopyNumber` and the predictor is
#' one specified feature column from the input data frame. It performs basic
#' cleaning by dropping rows with missing or non-finite values, and skips
#' fitting if too few usable observations remain.
#'
#' The function is intended for exploratory screening across multiple
#' candidate covariates to assess their linear association with `CopyNumber`.
#' Each call handles one predictor at a time, and returns either a fitted
#' `lm` object or `NULL` if fitting is not possible. The number of rows used
#' in the model is attached as an attribute (`n_used`) to the returned model
#' for downstream auditing.
#'
#' @param df A `data.frame` or tibble containing at least two numeric columns:
#'   \itemize{
#'     \item `CopyNumber`: the response variable (numeric).
#'     \item a predictor column specified by `xvar`.
#'   }
#' @param xvar A character scalar giving the name of the predictor variable
#'   to regress `CopyNumber` on.
#' @param min_n Integer (default = 20). Minimum number of complete rows required
#'   to attempt model fitting. If fewer rows are available after filtering,
#'   the function returns `NULL` with a warning.
#'
#' @return A fitted `lm` object (linear model) with formula
#'   `CopyNumber ~ xvar`, or `NULL` if the fit could not be performed.
#'   The fitted object includes an additional attribute, `"n_used"`, storing
#'   the number of observations used in the fit.
#'
#' @details
#' Steps performed by this function:
#' \enumerate{
#'   \item Selects only `CopyNumber` and the specified predictor column.
#'   \item Removes rows where either variable is missing or non-finite (Inf/NaN).
#'   \item If fewer than `min_n` rows remain, skips fitting.
#'   \item Attempts to fit a simple linear model using `lm()`. If fitting fails,
#'         returns `NULL` with a warning.
#' }
#'
#' This design is useful in pipelines where multiple features are screened
#' programmatically and robust handling of problematic variables is required.
#'
#' @examples
#' # Example data
#' df <- data.frame(
#'   CopyNumber = rnorm(100, mean = 5),
#'   TissueArea_mm = runif(100, 1, 10),
#'   pos_per_mm = rpois(100, 5)
#' )
#'
#' # Fit CopyNumber ~ TissueArea_mm
#' fit <- fit_lm_by_feature(df, "TissueArea_mm")
#' summary(fit)
#'
#' # Fit CopyNumber ~ pos_per_mm
#' fit2 <- fit_lm_by_feature(df, "pos_per_mm")
#' attr(fit2, "n_used")  # number of rows used
#'
#' @seealso [lm()], [audit()] for quickly checking completeness across multiple variables.
#'
#' @export
fit_lm_by_feature <- function(df, xvar, min_n = 20) {
  # keep only rows with finite CopyNumber and finite xvar
  df2 <- df %>%
    select(CopyNumber, all_of(xvar)) %>%
    filter(is.finite(CopyNumber), is.finite(.data[[xvar]])) %>%
    na.omit()
  
  n_used <- nrow(df2)
  if (n_used < min_n) {
    warning(sprintf("Skipping '%s': only %d usable rows after cleaning.", xvar, n_used))
    return(NULL)
  }
  
  # fit LM; if it still fails, return NULL
  out <- try(lm(as.formula(paste0("CopyNumber ~ ", xvar)), data = df2), silent = TRUE)
  if (inherits(out, "try-error")) {
    warning(sprintf("LM failed for '%s': %s", xvar, as.character(out)))
    return(NULL)
  }
  
  attr(out, "n_used") <- n_used
  out
}

# --- OPTIONAL: quick audit of finiteness before fitting ---
audit <- function(df, vars) {
  tibble(
    var = vars,
    n_nonfinite = sapply(vars, function(v)
      sum(!is.finite(df[[v]]), na.rm = TRUE) + sum(!is.finite(df$CopyNumber), na.rm = TRUE)),
    n_complete = sapply(vars, function(v)
      sum(is.finite(df[[v]]) & is.finite(df$CopyNumber), na.rm = TRUE))
  )
}


#' Compute Classification Performance Metrics from Confusion Matrix Counts
#' 
#' #' Customized for Estes SIV/HIV rebound study JDE006
#'
#' This function calculates a standard set of performance metrics given
#' the raw counts from a binary classification confusion matrix:
#' True Positives (TP), False Positives (FP), True Negatives (TN),
#' and False Negatives (FN). It is designed to be robust to small sample
#' sizes, missing denominators, and potential integer overflow.
#'
#' @param TP Numeric scalar. Count of true positives (predicted positive, actually positive).
#' @param FP Numeric scalar. Count of false positives (predicted positive, actually negative).
#' @param TN Numeric scalar. Count of true negatives (predicted negative, actually negative).
#' @param FN Numeric scalar. Count of false negatives (predicted negative, actually positive).
#'
#' @return A tibble with one row containing the following metrics:
#' \itemize{
#'   \item \strong{Sensitivity} (Recall, True Positive Rate): TP / (TP + FN).
#'   \item \strong{Specificity} (True Negative Rate): TN / (TN + FP).
#'   \item \strong{Precision} (Positive Predictive Value, PPV): TP / (TP + FP).
#'   \item \strong{NPV} (Negative Predictive Value): TN / (TN + FN).
#'   \item \strong{F1} score: harmonic mean of Precision and Sensitivity.
#'   \item \strong{Accuracy}: (TP + TN) / (TP + TN + FP + FN).
#'   \item \strong{BalancedAccuracy}: mean of Sensitivity and Specificity.
#'   \item \strong{YoudenJ}: Sensitivity + Specificity – 1, a threshold-free diagnostic index.
#'   \item \strong{MCC} (Matthews Correlation Coefficient): correlation-like statistic
#'         robust to imbalanced classes, with safeguards against overflow/undefined denominators.
#' }
#'
#' @details
#' All metrics are computed with checks to avoid division by zero or invalid
#' results. If a denominator is zero, the corresponding metric is returned as
#' `NA_real_`. The MCC denominator is protected against overflow by computing
#' on double precision and requiring finiteness and positivity before division.
#'
#' @examples
#' # Example confusion matrix counts
#' TP <- 69; FP <- 59; TN <- 484; FN <- 132
#' compute_metrics(TP, FP, TN, FN)
#'
#' # Edge case: zero positives
#' compute_metrics(TP = 0, FP = 10, TN = 90, FN = 0)
#'
#' @seealso [caret::confusionMatrix()] for generating confusion matrices directly.
#'
#' @export
compute_metrics <- function(TP, FP, TN, FN) {
  # cast to double early to avoid integer overflow
  TP <- as.numeric(TP); FP <- as.numeric(FP)
  TN <- as.numeric(TN); FN <- as.numeric(FN)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  acc  <- (TP + TN) / (TP + TN + FP + FN)
  
  f1 <- if (is.na(ppv) || is.na(sens) || (ppv + sens) == 0) NA_real_ else {
    2 * ppv * sens / (ppv + sens)
  }
  
  bacc   <- mean(c(sens, spec), na.rm = TRUE)
  youden <- sens + spec - 1
  
  # MCC with overflow-safe denominator
  denom <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
  denom <- if (is.finite(denom) && denom > 0) sqrt(denom) else NA_real_
  mcc   <- if (!is.na(denom) && denom > 0) ((TP * TN) - (FP * FN)) / denom else NA_real_
  
  tibble::tibble(
    Sensitivity = sens,
    Specificity = spec,
    Precision   = ppv,
    NPV         = npv,
    F1          = f1,
    Accuracy    = acc,
    BalancedAccuracy = bacc,
    YoudenJ     = youden,
    MCC         = mcc
  )
}