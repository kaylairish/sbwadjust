# Stabilized IPTW weights, moved into the package alongside boot_km_ratio()
# (CODE_PACKAGE_PLAN.md Step 4): boot_km_ratio()'s weight_type = "IPW" path
# calls this function directly, so once boot_km_ratio() lives in the package
# namespace it can no longer resolve a copy left behind in
# survival-risk-ratio/simulation_functions_survival.R (package functions
# don't search the caller's global environment). Previously only used by
# survival-risk-ratio; not part of Step 0's duplicate audit since it wasn't
# duplicated anywhere, just a dependency that has to move with its caller.

#' Stabilized inverse-probability-of-treatment weights for KM estimation
#'
#' `w = pi / ps` for treated units, `(1 - pi) / (1 - ps)` for control units,
#' where `ps` is the estimated propensity score and `pi` is the marginal
#' treatment probability used for stabilization.
#'
#' @param X_subset Covariate data frame for the full study sample.
#' @param A Treatment indicator (0/1).
#' @param use_glm If `TRUE` (default), estimate the propensity score via
#'   `stats::glm(A ~ ., family = binomial())`. If `FALSE`, use the fixed
#'   `p_known` for every unit (e.g. for a known randomization probability).
#' @param p_known Propensity score used for every unit when `use_glm = FALSE`.
#' @param pi Marginal treatment probability used for stabilization; if `NULL`
#'   (default), uses `mean(A)`.
#' @param trim Currently unused; reserved for propensity-score trimming.
#' @param return_ps Currently unused; reserved for returning the fitted
#'   propensity scores alongside the weights.
#' @return Numeric vector of stabilized IPTW weights, one per row of
#'   `X_subset`.
#' @export
get_ipws_for_study = function(X_subset, A,
                              use_glm = TRUE,
                              p_known = 0.5,           # used only when use_glm = FALSE
                              pi = NULL,               # if NULL, uses mean(A)
                              trim = c(0.01, 0.99),
                              return_ps = FALSE) {
  A = as.integer(A)

  # propensity score
  if (!use_glm) {
    ps = rep(p_known, length(A))
  } else {
    df = as.data.frame(X_subset)
    df$A = A
    fit = stats::glm(A ~ ., family = stats::binomial(), data = df)
    ps = as.numeric(stats::predict(fit, type = "response"))
  }

  # marginal treatment probability for stabilization
  if (is.null(pi)) pi = mean(A)

  # stabilized IPTW
  w = ifelse(A == 1, pi / ps, (1 - pi) / (1 - ps))

  w
}
