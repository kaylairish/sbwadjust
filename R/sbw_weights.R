# User-facing formula API (design-stage): sbw_weights() wraps
# get_sbws_for_study() (weights.R) behind a formula/data/treatment interface,
# with S3 methods (print, summary, plot, weights) on the fitted object.
# Added for the v0.2 "usable by a stranger" release described in
# manuscript/sbwadjust_package_design_note.tex.

#' Build the balance design matrix from a one-sided formula
#'
#' @param balance One-sided formula naming the covariates to balance.
#' @param data Data frame containing those covariates.
#' @return A data frame of numeric balance columns (factors dummy-coded,
#'   intercept dropped).
#' @keywords internal
.build_design = function(balance, data) {
  if (!inherits(balance, "formula") || length(balance) != 2L) {
    stop("`balance` must be a one-sided formula, e.g. ~ age + sex + bmi.")
  }
  mf = stats::model.frame(balance, data = data, na.action = stats::na.pass)
  mm = stats::model.matrix(balance, data = mf)
  keep = colnames(mm) != "(Intercept)"
  as.data.frame(mm[, keep, drop = FALSE])
}

#' Recode a raw treatment vector to 0/1, recording the level map used
#'
#' @param raw Raw treatment vector: numeric 0/1, logical, or a two-level
#'   factor/character vector.
#' @return A list with `A` (integer 0/1 vector) and `levels` (named integer
#'   vector mapping each original level/label to 0 or 1; for numeric input the
#'   names are "0" and "1", and for logical input "FALSE" and "TRUE").
#' @keywords internal
.recode_treatment = function(raw) {
  if (is.logical(raw)) {
    return(list(A = as.integer(raw), levels = c(`FALSE` = 0L, `TRUE` = 1L)))
  }
  if (is.numeric(raw)) {
    u = sort(unique(raw))
    if (!isTRUE(all.equal(u, c(0, 1)))) {
      stop("Numeric `treatment` must be coded 0/1.")
    }
    return(list(A = as.integer(raw), levels = c(`0` = 0L, `1` = 1L)))
  }
  f = factor(raw)
  lv = levels(f)
  if (length(lv) != 2L) stop("`treatment` must have exactly two levels.")
  levels_map = stats::setNames(c(0L, 1L), lv)
  list(A = as.integer(levels_map[as.character(f)]), levels = levels_map)
}

#' Resolve the `treatment` argument (bare column name or vector) against `data`
#'
#' @param treatment_expr The unevaluated `treatment` argument, from
#'   `substitute(treatment)`.
#' @param data Data frame to resolve a bare column name against.
#' @param env Environment to evaluate `treatment_expr` in if it isn't a
#'   column of `data` (e.g. a vector supplied directly).
#' @return A list with `name` (character, for reporting/bootstrap re-lookup)
#'   and `raw` (the resolved vector).
#' @keywords internal
.resolve_treatment = function(treatment_expr, data, env) {
  if (is.symbol(treatment_expr) && as.character(treatment_expr) %in% names(data)) {
    name = as.character(treatment_expr)
    return(list(name = name, raw = data[[name]]))
  }
  raw = eval(treatment_expr, data, env)
  name = if (is.symbol(treatment_expr)) as.character(treatment_expr) else deparse(treatment_expr)
  list(name = name, raw = raw)
}

#' Fit stable balancing weights for a two-arm study (design stage)
#'
#' Formula-based front end to [get_sbws_for_study()]: balances the covariates
#' named in `balance` (factors dummy-coded automatically) to the pooled-arm
#' mean, exactly (imbalance tolerance fixed at zero), matching the paper's
#' SBW specification. Intended to be called before unblinding, using only
#' baseline covariates and the randomization assignment.
#'
#' @param balance One-sided formula naming the covariates to balance, e.g.
#'   `~ age + sex + bmi + region`.
#' @param data Data frame containing the balance covariates and treatment.
#' @param treatment Treatment assignment: either a bare column name in
#'   `data` (unquoted, as in `treatment = arm`) or a vector. Accepts numeric
#'   0/1, logical, or a two-level factor/character vector.
#' @return An object of class `sbw_fit` with elements `weights`, `n_clipped`,
#'   `max_abs_clipped` (see [get_sbws_for_study()]), `balance` (the formula),
#'   `treatment_name`, `treatment_levels`, `treatment` (recoded 0/1),
#'   `X` (the balance design matrix), `data`, `n`, and `call`.
#' @examples
#' set.seed(1)
#' n = 100
#' trial_baseline = data.frame(
#'   age = rnorm(n, 50, 10),
#'   region = sample(c("N", "S"), n, replace = TRUE),
#'   arm = rbinom(n, 1, 0.5)
#' )
#' sbw = sbw_weights(~ age + region, data = trial_baseline, treatment = arm)
#' sbw
#' summary(sbw)
#' @export
sbw_weights = function(balance, data, treatment) {
  data = as.data.frame(data)
  treatment_expr = substitute(treatment)
  resolved = .resolve_treatment(treatment_expr, data, parent.frame())
  rec = .recode_treatment(resolved$raw)
  A = rec$A

  X = .build_design(balance, data)
  if (anyNA(X)) stop("Missing values in the balance covariates are not supported.")
  if (nrow(X) != length(A)) stop("`data` and `treatment` do not have the same length.")

  fit = get_sbws_for_study(X, A)

  structure(
    list(
      weights = fit$w,
      n_clipped = fit$n_clipped,
      max_abs_clipped = fit$max_abs_clipped,
      balance = balance,
      treatment_name = resolved$name,
      treatment_levels = rec$levels,
      treatment = A,
      X = X,
      data = data,
      n = nrow(data),
      call = match.call()
    ),
    class = "sbw_fit"
  )
}

#' @export
print.sbw_fit = function(x, ...) {
  cat("<sbw_fit>\n")
  cat("  balance: ", deparse(x$balance), "\n", sep = "")
  cat("  n:       ", x$n, " (", sum(x$treatment == 1L), " treated, ",
      sum(x$treatment == 0L), " control)\n", sep = "")
  if (x$n_clipped > 0L) {
    cat("  note:    ", x$n_clipped,
        " weight(s) clipped at 0 (nonnegative-QP fallback used)\n", sep = "")
  }
  invisible(x)
}

#' @importFrom stats weights
#' @export
weights.sbw_fit = function(object, ...) object$weights

#' Summarize an `sbw_fit`: balance table and weight diagnostics
#'
#' @param object An `sbw_fit` from [sbw_weights()].
#' @param ... Currently unused.
#' @return An object of class `summary.sbw_fit`, printed by
#'   `print.summary.sbw_fit()`.
#' @export
summary.sbw_fit = function(object, ...) {
  X = object$X
  A = object$treatment
  w = object$weights

  wmean = function(v, wt) sum(v * wt) / sum(wt)
  ess = function(wt) sum(wt)^2 / sum(wt^2)

  balance_tbl = data.frame(
    covariate = colnames(X),
    mean_treated_unw = vapply(X, function(col) mean(col[A == 1L]), numeric(1)),
    mean_control_unw = vapply(X, function(col) mean(col[A == 0L]), numeric(1)),
    mean_treated_w = vapply(X, function(col) wmean(col[A == 1L], w[A == 1L]), numeric(1)),
    mean_control_w = vapply(X, function(col) wmean(col[A == 0L], w[A == 0L]), numeric(1)),
    row.names = NULL
  )

  structure(
    list(
      balance = balance_tbl,
      ess_treated = ess(w[A == 1L]),
      ess_control = ess(w[A == 0L]),
      n_treated = sum(A == 1L),
      n_control = sum(A == 0L),
      n_clipped = object$n_clipped,
      max_abs_clipped = object$max_abs_clipped
    ),
    class = "summary.sbw_fit"
  )
}

#' @export
print.summary.sbw_fit = function(x, ...) {
  cat("Balance (unweighted vs. SBW-weighted arm means):\n")
  print(x$balance, row.names = FALSE)
  cat("\nEffective sample size: ",
      round(x$ess_treated, 1), " treated (of ", x$n_treated, "), ",
      round(x$ess_control, 1), " control (of ", x$n_control, ")\n", sep = "")
  if (x$n_clipped > 0L) {
    cat("Clipped weights: ", x$n_clipped,
        " (max magnitude ", signif(x$max_abs_clipped, 3), ")\n", sep = "")
  }
  invisible(x)
}

#' Plot covariate balance before and after SBW weighting
#'
#' Standardized mean differences (treated vs. control) for each balance
#' covariate, unweighted and SBW-weighted.
#'
#' @param x An `sbw_fit` from [sbw_weights()].
#' @param ... Passed on to the underlying `graphics::plot()` call.
#' @return `x`, invisibly.
#' @export
plot.sbw_fit = function(x, ...) {
  X = x$X
  A = x$treatment
  w = x$weights

  wmean = function(v, wt) sum(v * wt) / sum(wt)
  wvar = function(v, wt) sum(wt * (v - wmean(v, wt))^2) / sum(wt)
  smd = function(v1, v0) (mean(v1) - mean(v0)) / sqrt((stats::var(v1) + stats::var(v0)) / 2)
  wsmd = function(v1, wt1, v0, wt0) {
    (wmean(v1, wt1) - wmean(v0, wt0)) / sqrt((wvar(v1, wt1) + wvar(v0, wt0)) / 2)
  }

  unw = vapply(X, function(col) smd(col[A == 1L], col[A == 0L]), numeric(1))
  wtd = vapply(X, function(col) {
    wsmd(col[A == 1L], w[A == 1L], col[A == 0L], w[A == 0L])
  }, numeric(1))

  covs = colnames(X)
  yy = seq_along(covs)
  xr = range(c(unw, wtd, 0), na.rm = TRUE)

  graphics::plot(unw, yy, pch = 1, xlim = xr, yaxt = "n",
                 xlab = "Standardized mean difference", ylab = "",
                 main = "Covariate balance", ...)
  graphics::points(wtd, yy, pch = 16)
  graphics::axis(2, at = yy, labels = covs, las = 1)
  graphics::abline(v = 0, lty = 2)
  graphics::legend("topright", legend = c("Unweighted", "SBW-weighted"),
                    pch = c(1, 16), bty = "n")
  invisible(x)
}
