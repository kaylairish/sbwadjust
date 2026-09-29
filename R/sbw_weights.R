# User-facing formula API (design-stage): sbw_weights() wraps
# get_sbws_for_study() (weights.R) behind a formula/data/treatment interface,
# with S3 methods (print, summary, plot, weights) on the fitted object.

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
  if (anyNA(raw)) stop("Missing values in `treatment` are not supported.")
  if (length(unique(raw)) < 2L) {
    stop("`treatment` has only one arm; both arms need at least one unit.")
  }
  if (is.logical(raw)) {
    return(list(A = as.integer(raw), levels = c(`FALSE` = 0L, `TRUE` = 1L)))
  }
  if (is.numeric(raw)) {
    if (!all(raw %in% c(0, 1))) stop("Numeric `treatment` must be coded 0/1.")
    return(list(A = as.integer(raw), levels = c(`0` = 0L, `1` = 1L)))
  }
  f = factor(raw)
  lv = levels(f)
  if (length(lv) != 2L) stop("`treatment` must have exactly two levels.")
  levels_map = stats::setNames(c(0L, 1L), lv)
  list(A = as.integer(levels_map[as.character(f)]), levels = levels_map)
}

#' Resolve the `treatment` argument (column name or vector) against `data`
#'
#' @param treatment_expr The unevaluated `treatment` argument, from
#'   `substitute(treatment)`.
#' @param data Data frame to resolve a column name against.
#' @param env Environment to evaluate `treatment_expr` in if it isn't a
#'   column of `data` (e.g. a vector supplied directly).
#' @return A list with `name` (character, a label for reporting) and `raw`
#'   (the resolved vector).
#' @keywords internal
.resolve_treatment = function(treatment_expr, data, env) {
  if (is.symbol(treatment_expr) && as.character(treatment_expr) %in% names(data)) {
    name = as.character(treatment_expr)
    return(list(name = name, raw = data[[name]]))
  }
  raw = eval(treatment_expr, data, env)
  # a single string is a column name, quoted (treatment = "arm") or held in a
  # variable (treatment = trt_col), not a one-element treatment vector
  if (is.character(raw) && length(raw) == 1L) {
    if (!raw %in% names(data)) {
      stop("`treatment` \"", raw, "\" is not a column of `data`.")
    }
    return(list(name = raw, raw = data[[raw]]))
  }
  name = if (is.symbol(treatment_expr)) {
    as.character(treatment_expr)
  } else {
    paste(deparse(treatment_expr, width.cutoff = 500L), collapse = " ")
  }
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
#' @param treatment Treatment assignment: either a column name in `data`,
#'   unquoted (`treatment = arm`) or quoted (`treatment = "arm"`), or a vector. Accepts numeric
#'   0/1, logical, or a two-level factor/character vector. For a factor, the
#'   second level is treated (coded 1); for a character vector, the second in
#'   alphabetical order. Use a factor with explicit `levels` to control this.
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

  # the closed-form solve needs each arm's [1, X] to have full column rank
  for (arm in c(1L, 0L)) {
    Xa = cbind(1, as.matrix(X[A == arm, , drop = FALSE]))
    if (qr(Xa)$rank < ncol(Xa)) {
      stop("Balance covariates are collinear within the ",
           if (arm == 1L) "treated" else "control",
           " arm; drop or combine redundant ones.")
    }
  }

  # translate solver errors that the rank check doesn't catch
  fit = tryCatch(get_sbws_for_study(X, A), error = function(e) e)
  if (inherits(fit, "error")) {
    msg = conditionMessage(fit)
    if (grepl("constraints are inconsistent", msg, fixed = TRUE)) {
      stop("Exact balance is infeasible: no nonnegative weights in one arm ",
           "reproduce the pooled mean of the balance covariates. This is an ",
           "overlap problem: on some covariate (or combination of covariates), ",
           "one arm's values all lie on one side of the pooled mean. Compare ",
           "the arms' ranges, e.g. tapply(x, arm, range).")
    }
    if (grepl("singular", msg, fixed = TRUE)) {
      stop("Balance covariates are collinear within one arm; drop or combine redundant ones.")
    }
    stop(fit)
  }

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
  cat("  balance: ", paste(deparse(x$balance, width.cutoff = 500L), collapse = " "), "\n", sep = "")
  # name the arms when the treated level isn't self-evident (factor/character)
  lv = names(x$treatment_levels)
  arm_label = if (identical(lv, c("0", "1")) || identical(lv, c("FALSE", "TRUE"))) {
    c("", "")
  } else {
    paste0(" [", lv, "]")
  }
  cat("  n:       ", x$n, " (", sum(x$treatment == 1L), " treated", arm_label[2], ", ",
      sum(x$treatment == 0L), " control", arm_label[1], ")\n", sep = "")
  # count zero weights directly: n_clipped counts QP rounding noise, not units
  # dropped, and varies by platform
  zero = x$weights < 1e-10
  n_arm = c(treated = sum(x$treatment == 1L), control = sum(x$treatment == 0L))
  n_zero = c(treated = sum(zero & x$treatment == 1L), control = sum(zero & x$treatment == 0L))
  if (any(n_zero > 0L)) {
    parts = paste0(n_zero, " of ", n_arm, " ", names(n_arm))[n_zero > 0L]
    cat("  note:    ", paste(parts, collapse = " and "),
        " units got weight 0, so they don't contribute\n",
        "           to the estimate; the rest carry the balance. See summary().\n", sep = "")
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
