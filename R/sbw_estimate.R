# User-facing analysis-stage API: sbw_estimate() dispatches to one of a
# closed menu of estimands (design note's "no user-supplied functional"
# principle -- an uncovered estimand means take weights(sbw) and run your
# own estimator). ATE/RR/mann_whitney/quantile_diff/quantile_ratio share one
# bootstrap engine here (resample rows, refit SBW weights, recompute the
# point estimate); survival_ratio instead calls the already-tested
# boot_km_ratio() from km_ratio.R rather than duplicating it.
#
# Estimand coverage note (per manuscript/main.tex and
# main_jasa_supplement_body.tex): ATE, RR, and survival_ratio have full
# worked estimators/simulations in the paper. mann_whitney implements only
# the finite/uncensored case given in main_jasa_supplement_body.tex (the
# right-censored case is left to future work, per that section). RMST is
# intentionally not included in this release. quantile_diff/quantile_ratio
# (main.tex's "Quantiles and medians" item) use the natural SBW-weighted
# empirical-quantile plug-in -- the generalized-inverse weighted quantile,
# arm-specific difference/ratio, bootstrap CI -- which is not spelled out
# verbatim in the paper as an estimator recipe; flagged for review.

#' Weighted mean
#' @keywords internal
.weighted_mean = function(x, w) sum(x * w) / sum(w)

#' Weighted generalized-inverse (step-function) quantile
#'
#' `q_hat(u) = inf{y : F_w(y) >= u}` for the weighted empirical CDF `F_w`,
#' the natural SBW-weighted analogue of `stats::quantile(type = 1)`.
#'
#' Normalizes by `sum(w)`, not `length(y)`. By construction (the balancing
#' QP's intercept equality constraint), SBW weights for an arm sum to that
#' arm's size -- but only up to solver precision, since the nonneg-QP
#' fallback's clipping of negative numerical dust can nudge the sum a hair
#' away from it. Dividing by `sum(w)` keeps `F_w` a valid CDF (`F_w(Inf) = 1`
#' exactly) regardless, and avoids having to track which arm's size applies
#' to whatever `w` subset was passed in.
#'
#' @param y Numeric outcome vector for one arm.
#' @param w Nonnegative weights, same length as `y`.
#' @param probs Vector of probabilities in (0, 1).
#' @return Numeric vector of estimated quantiles, one per element of `probs`.
#' @keywords internal
.weighted_quantile = function(y, w, probs) {
  ord = order(y)
  y_sorted = y[ord]
  cw = cumsum(w[ord]) / sum(w)
  vapply(probs, function(u) y_sorted[which(cw >= u)[1]], numeric(1))
}

#' Weighted Mann-Whitney / win-probability estimator
#'
#' `P(Y1 > Y0) + 0.5 P(Y1 = Y0)`, weighted by SBW weights within each arm.
#' Finite/uncensored outcomes only.
#'
#' @param y1,w1 Outcome and weights for the treated arm.
#' @param y0,w0 Outcome and weights for the control arm.
#' @return Scalar win-probability estimate.
#' @keywords internal
.weighted_win_prob = function(y1, w1, y0, w0) {
  cmp = outer(y1, y0, function(a, b) (a > b) + 0.5 * (a == b))
  sum(cmp * outer(w1, w0)) / (sum(w1) * sum(w0))
}

#' Estimate a treatment effect from a fitted `sbw_fit` object
#'
#' @param object An `sbw_fit` from [sbw_weights()].
#' @param outcome A one-sided-response formula naming the outcome, e.g.
#'   `Y ~ 1` (or `Surv(time, status) ~ 1` for `estimand = "survival_ratio"`,
#'   which requires the `survival` package to be attached so `Surv()`
#'   resolves in the formula's environment).
#' @param estimand One of `"ATE"`, `"RR"` (relative risk, i.e. ratio of
#'   weighted arm means), `"survival_ratio"`, `"mann_whitney"`,
#'   `"quantile_diff"`, or `"quantile_ratio"`. No user-supplied functional is
#'   accepted; an uncovered estimand means: take `weights(object)` and run
#'   (and bootstrap) your own estimator.
#' @param ... Passed to methods.
#' @return An object of class `sbw_estimate`.
#' @export
sbw_estimate = function(object, outcome, estimand, ...) {
  UseMethod("sbw_estimate")
}

#' @param probs Vector of probabilities in (0, 1), used when `estimand` is
#'   `"quantile_diff"` or `"quantile_ratio"`. Default `0.5` (the median).
#' @param horizon Time `t0` at which to evaluate the survival ratio; required
#'   for `estimand = "survival_ratio"`.
#' @param B Number of bootstrap replicates.
#' @param alpha Significance level for the confidence interval.
#' @param ci_method Either `"wald"` (default) or `"percentile"`.
#' @param seed Optional seed set at the start of the bootstrap.
#' @examples
#' set.seed(1)
#' n = 100
#' trial_baseline = data.frame(
#'   age = rnorm(n, 50, 10),
#'   region = sample(c("N", "S"), n, replace = TRUE),
#'   arm = rbinom(n, 1, 0.5)
#' )
#' trial_baseline$Y = rbinom(n, 1, plogis(-1 + 0.02 * trial_baseline$age))
#' sbw = sbw_weights(~ age + region, data = trial_baseline, treatment = arm)
#' sbw_estimate(sbw, Y ~ 1, estimand = "RR", B = 200, seed = 1)
#' @rdname sbw_estimate
#' @export
sbw_estimate.sbw_fit = function(object, outcome, estimand,
                                 probs = 0.5, horizon = NULL,
                                 B = 1500, alpha = 0.05,
                                 ci_method = c("wald", "percentile"),
                                 seed = NULL, ...) {
  estimand = match.arg(estimand, c("ATE", "RR", "survival_ratio",
                                    "mann_whitney", "quantile_diff", "quantile_ratio"))
  ci_method = match.arg(ci_method)
  data = object$data
  A = object$treatment

  if (estimand == "survival_ratio") {
    if (is.null(horizon)) {
      stop("`horizon` (time t0) is required for estimand = \"survival_ratio\".")
    }
    mf = stats::model.frame(outcome, data = data)
    resp = stats::model.response(mf)
    if (!inherits(resp, "Surv")) {
      stop("`outcome` must be a Surv(time, status) formula for estimand = \"survival_ratio\".")
    }
    time = resp[, 1]
    status = resp[, 2]

    res = boot_km_ratio(time, status, A, object$X, t0 = horizon,
                         B = B, alpha = alpha, weight_type = "SBW",
                         ci_method = ci_method, seed = seed)

    return(structure(
      list(
        estimand = "survival_ratio", scale = "log",
        estimate = stats::setNames(res$est, "survival_ratio"),
        se = res$se_log, ci = matrix(res$ci, ncol = 2),
        alpha = alpha, B = B, horizon = horizon,
        mc_fail = res$MC_fail, detail = res
      ),
      class = "sbw_estimate"
    ))
  }

  mf = stats::model.frame(outcome, data = data, na.action = stats::na.pass)
  Y = stats::model.response(mf)
  if (anyNA(Y)) stop("Missing values in the outcome are not supported.")

  point_fun = switch(estimand,
    ATE = function(w, A, Y) {
      .weighted_mean(Y[A == 1L], w[A == 1L]) - .weighted_mean(Y[A == 0L], w[A == 0L])
    },
    RR = function(w, A, Y) {
      log(.weighted_mean(Y[A == 1L], w[A == 1L]) / .weighted_mean(Y[A == 0L], w[A == 0L]))
    },
    mann_whitney = function(w, A, Y) {
      .weighted_win_prob(Y[A == 1L], w[A == 1L], Y[A == 0L], w[A == 0L])
    },
    quantile_diff = function(w, A, Y) {
      .weighted_quantile(Y[A == 1L], w[A == 1L], probs) -
        .weighted_quantile(Y[A == 0L], w[A == 0L], probs)
    },
    quantile_ratio = function(w, A, Y) {
      log(.weighted_quantile(Y[A == 1L], w[A == 1L], probs) /
            .weighted_quantile(Y[A == 0L], w[A == 0L], probs))
    }
  )

  est = point_fun(object$weights, A, Y)

  n = object$n
  if (!is.null(seed)) set.seed(seed)
  boot_reps = matrix(NA_real_, nrow = B, ncol = length(est))
  for (b in seq_len(B)) {
    idx = sample.int(n, n, replace = TRUE)
    data_b = data[idx, , drop = FALSE]
    rep_est = tryCatch({
      X_b = .build_design(object$balance, data_b)
      A_b = .apply_treatment_levels(data_b[[object$treatment_name]], object$treatment_levels)
      w_b = get_sbws_for_study(X_b, A_b)$w
      mf_b = stats::model.frame(outcome, data = data_b, na.action = stats::na.pass)
      Y_b = stats::model.response(mf_b)
      point_fun(w_b, A_b, Y_b)
    }, error = function(e) rep(NA_real_, length(est)))
    boot_reps[b, ] = rep_est
  }

  finite_rows = stats::complete.cases(boot_reps)
  if (sum(finite_rows) < 2L) {
    stop("Bootstrap failed for nearly all resamples; cannot estimate a standard error.")
  }
  se = apply(boot_reps[finite_rows, , drop = FALSE], 2, stats::sd)
  z = stats::qnorm(1 - alpha / 2)

  if (ci_method == "wald") {
    ci = cbind(est - z * se, est + z * se)
  } else {
    ci = t(apply(boot_reps[finite_rows, , drop = FALSE], 2, stats::quantile,
                 probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  }

  scale = if (estimand %in% c("RR", "quantile_ratio")) "log" else "identity"
  estimate_out = if (scale == "log") exp(est) else est
  ci_out = if (scale == "log") exp(ci) else ci
  ci_out = matrix(ci_out, ncol = 2)

  nm = if (estimand %in% c("quantile_diff", "quantile_ratio")) paste0("q", probs) else estimand
  names(estimate_out) = nm

  structure(
    list(
      estimand = estimand, scale = scale,
      estimate = estimate_out, se = se, ci = ci_out,
      alpha = alpha, B = B,
      boot_fail_rate = mean(!finite_rows)
    ),
    class = "sbw_estimate"
  )
}

#' @export
print.sbw_estimate = function(x, ...) {
  cat("<sbw_estimate: ", x$estimand, ">\n", sep = "")
  est = x$estimate
  ci = x$ci
  nm = names(est)
  if (is.null(nm)) nm = rep(x$estimand, length(est))
  for (i in seq_along(est)) {
    cat(sprintf("  %s: %.4g  (%.0f%% CI: %.4g, %.4g)\n",
                nm[i], est[i], 100 * (1 - x$alpha), ci[i, 1], ci[i, 2]))
  }
  invisible(x)
}
