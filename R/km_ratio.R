# Shared KM-ratio estimator: S1(t0)/S0(t0) from (weighted) KM + log-log
# Greenwood CI, plus its bootstrap wrapper.
#
# Canonical picks (2026-07-22, CODE_PACKAGE_PLAN.md Step 2, confirmed with the author):
# - km_ratio_loglog_greenwood: no eps-clamping of S0/S1 before the log-log
#   transform. R's own `survival` package treats S=0/1 as genuinely undefined
#   for the log-log CI (returns NA there) rather than picking an arbitrary
#   epsilon, and this codebase's own saved results (checked across all
#   ~1M rows in survival-risk-ratio/results-power + results-t1e) never
#   actually hit that boundary, so dropping the clamp changes nothing already
#   reported.
# - boot_km_ratio: survival-risk-ratio's fuller version (IPW option, hard-fail-
#   to-unadjusted fallback with MC_fail, SBW clipping diagnostics), with an
#   added ci_method switch so illustration-of-approach/table_1.R keeps its
#   percentile CI instead of switching to Wald.

#' Extract survival + SE at a fixed time from a `survfit` object, by arm
#'
#' @param fit A `survival::survfit` object fit on `Surv(time, status) ~ A`.
#' @param t0 Time at which to extract survival estimates.
#' @return A list with `S0`, `S1` (survival estimates for A=0/A=1) and `se0`,
#'   `se1` (their standard errors).
#' @keywords internal
.extract_surv_at_t0 = function(fit, t0) {
  s = summary(fit, times = t0, extend = TRUE)
  strata = as.character(s$strata)
  S      = as.numeric(s$surv)
  seS    = as.numeric(s$std.err)

  find_arm = function(val) {
    idx = grep(paste0("A=", val, "\\b"), strata)
    if (length(idx) == 0) idx = grep(paste0("A=", val), strata)
    if (length(idx) == 0) stop("Couldn't identify KM strata for A=", val)
    idx[1]
  }

  i0 = find_arm(0)
  i1 = find_arm(1)

  list(S0 = S[i0], S1 = S[i1], se0 = seS[i0], se1 = seS[i1])
}

#' Weighted KM survival-ratio S1(t0)/S0(t0) with a log-log Greenwood CI
#'
#' @param time Event/censoring time.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param A Treatment indicator (0/1).
#' @param t0 Time at which to evaluate the survival ratio.
#' @param alpha Significance level for the confidence interval (default 0.05).
#' @param weights Optional case weights (e.g. SBW or IPW), passed to
#'   `survival::survfit()`.
#' @return A list with `S0`, `S1`, `ratio` (= S1(t0)/S0(t0)), `log_ratio`,
#'   `se_log`, `ci_ratio`, and `ci_log`.
#' @export
km_ratio_loglog_greenwood = function(time, status, A, t0, alpha = 0.05, weights = NULL) {
  df = data.frame(time = time, status = status, A = as.integer(A))

  fit = survival::survfit(survival::Surv(time, status) ~ A, data = df, weights = weights)

  out = .extract_surv_at_t0(fit, t0 = t0)
  S0 = out$S0
  S1 = out$S1
  se0 = out$se0
  se1 = out$se1

  # log-log Greenwood delta
  gprime = function(S) 1 / (S * log(S))
  var_g0 = (gprime(S0) * se0)^2
  var_g1 = (gprime(S1) * se1)^2

  var_logS0 = (log(S0))^2 * var_g0
  var_logS1 = (log(S1))^2 * var_g1

  logRR = log(S1) - log(S0)
  se_logRR = sqrt(var_logS0 + var_logS1)

  z = stats::qnorm(1 - alpha / 2)
  ci_log = c(logRR - z * se_logRR, logRR + z * se_logRR)

  list(
    S0 = S0, S1 = S1,
    ratio = exp(logRR),          # this is SRR(t0)
    log_ratio = logRR,
    se_log = se_logRR,
    ci_ratio = exp(ci_log),
    ci_log = ci_log
  )
}

#' Bootstrap CI for a weighted-KM survival ratio (SBW or IPW)
#'
#' Bootstraps the standard error of [km_ratio_loglog_greenwood()]'s
#' log-ratio, then builds either a Wald CI (log scale) or a percentile CI
#' (log scale). If the point estimate or bootstrap SE come out non-finite,
#' falls back to the unadjusted KM ratio and sets `MC_fail = TRUE`; if even
#' the unadjusted ratio is undefined (a double-degenerate case with no valid
#' fallback), throws an error rather than returning a nonsense finite
#' estimate.
#'
#' @param time Event/censoring time.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param A Treatment indicator (0/1).
#' @param X_subset Covariate data frame for the full study sample, used to
#'   fit the SBW/IPW weights.
#' @param t0 Time at which to evaluate the survival ratio.
#' @param B Number of bootstrap replicates.
#' @param alpha Significance level for the confidence interval (default 0.05).
#' @param weight_type Either `"SBW"` ([get_sbws_for_study()]) or `"IPW"`
#'   ([get_ipws_for_study()]).
#' @param ci_method Either `"wald"` (default) or `"percentile"`.
#' @param ipw_use_glm Passed through to [get_ipws_for_study()] when
#'   `weight_type = "IPW"`.
#' @param seed Optional seed set at the start of the bootstrap.
#' @param verbose Currently unused; reserved for future diagnostic output.
#' @return A list with `MC_fail`, `log_est`, `est`, `se_log`, `ci_log`, `ci`,
#'   `boot_fail_rate`, `boot_n_finite_reps`, and SBW clipping diagnostics
#'   (`sbw_n_clipped_full`, `sbw_max_abs_clipped_full`,
#'   `sbw_n_clipped_boot_mean`, `sbw_n_clipped_boot_max`,
#'   `sbw_max_abs_clipped_boot_max`; all `NA` when `weight_type = "IPW"`).
#' @export
boot_km_ratio = function(time, status, A, X_subset, t0,
                         B = 1500,
                         alpha = 0.05,
                         weight_type = c("SBW", "IPW"),
                         ci_method = c("wald", "percentile"),
                         ipw_use_glm = TRUE,
                         seed = NULL,
                         verbose = FALSE) {

  weight_type = match.arg(weight_type)
  ci_method = match.arg(ci_method)
  if (!is.null(seed)) set.seed(seed)

  n = length(A)
  z = stats::qnorm(1 - alpha/2)

  # -------------------------
  # 1) Point estimate attempt
  # -------------------------
  point_est = tryCatch({
    if (weight_type == "SBW") {
      sbw_full = get_sbws_for_study(X_subset, A)
      w_full = sbw_full$w
      nclip_full = sbw_full$n_clipped
      maxclip_full = sbw_full$max_abs_clipped
    } else {
      w_full = get_ipws_for_study(X_subset, A, use_glm = ipw_use_glm)
      nclip_full = NA_integer_
      maxclip_full = NA_real_
    }

    full = km_ratio_loglog_greenwood(time, status, A, t0 = t0, alpha = alpha, weights = w_full)

    list(
      finite_rep = is.finite(full$log_ratio),
      w_full = w_full,
      loghat = full$log_ratio,
      est = full$ratio,
      nclip_full = nclip_full,
      maxclip_full = maxclip_full
    )
  }, error = function(e) {
    list(finite_rep = FALSE)
  })

  if (!isTRUE(point_est$finite_rep)) {
    # HARD FAIL, fallback to unadjusted, mark MC_fail

    # unadjusted
    unadj = tryCatch(
      km_ratio_loglog_greenwood(time, status, A, t0 = t0, alpha = alpha, weights = NULL),
      error = function(e) NULL
    )
    if (is.null(unadj) || !is.finite(unadj$log_ratio)) {
      stop("Unadjusted KM ratio failed (cannot fallback).")
    }

    return(list(
      MC_fail = TRUE,
      log_est = unadj$log_ratio,
      est = unadj$ratio,
      se_log = unadj$se_log,
      ci_log = unadj$ci_log,
      ci = unadj$ci_ratio,
      boot_fail_rate = NA_real_,
      boot_n_finite_reps = NA_integer_,
      sbw_n_clipped_full = NA_integer_,
      sbw_max_abs_clipped_full = NA_real_,
      sbw_n_clipped_boot_mean = NA_real_,
      sbw_n_clipped_boot_max  = NA_integer_,
      sbw_max_abs_clipped_boot_max = NA_real_
    ))
  }

  loghat = point_est$loghat

  # -------------------------
  # 2) Bootstrap attempt
  # -------------------------
  boot_log = rep(NA_real_, B)
  fail = rep(FALSE, B)

  nclip_boot = if (weight_type == "SBW") rep(NA_integer_, B) else NULL
  maxclip_boot = if (weight_type == "SBW") rep(NA_real_, B) else NULL

  for (b in seq_len(B)) {
    idx = sample.int(n, n, replace = TRUE)
    time_boot = time[idx]
    status_boot = status[idx]
    A_boot = A[idx]
    X_boot = X_subset[idx, , drop = FALSE]

    bootstrap_point_est = tryCatch({
      if (weight_type == "SBW") {
        sbw_boot = get_sbws_for_study(X_boot, A_boot)
        w_boot = sbw_boot$w
        nclip_boot[b] = sbw_boot$n_clipped
        maxclip_boot[b] = sbw_boot$max_abs_clipped
      } else {
        w_boot = get_ipws_for_study(X_boot, A_boot, use_glm = ipw_use_glm)
      }

      km_ratio_loglog_greenwood(time_boot, status_boot, A_boot, t0 = t0, alpha = alpha, weights = w_boot)$log_ratio
    }, error = function(e) NA_real_)

    if (!is.finite(bootstrap_point_est)) {
      fail[b] = TRUE
      # use unadjusted for that bootstrap draw
      bootstrap_point_est = tryCatch(
        km_ratio_loglog_greenwood(time_boot, status_boot, A_boot, t0 = t0, alpha = alpha, weights = NULL)$log_ratio,
        error = function(e) NA_real_
      )
    }
    boot_log[b] = bootstrap_point_est
  }

  finite_rep = is.finite(boot_log)
  boot_n_finite_reps = sum(finite_rep)
  fail_rate = mean(fail)

  se_boot = stats::sd(boot_log[finite_rep])

  if (!is.finite(se_boot) || se_boot <= 0) {
    # unadjusted
    unadj = tryCatch(
      km_ratio_loglog_greenwood(time, status, A, t0 = t0, alpha = alpha, weights = NULL),
      error = function(e) NULL
    )
    if (is.null(unadj) || !is.finite(unadj$log_ratio)) {
      stop("Unadjusted KM ratio failed (cannot fallback).")
    }

    return(list(
      MC_fail = TRUE,
      log_est = unadj$log_ratio,
      est = unadj$ratio,
      se_log = unadj$se_log,
      ci_log = unadj$ci_log,
      ci = unadj$ci_ratio,
      boot_fail_rate = fail_rate,
      boot_n_finite_reps = boot_n_finite_reps,

      sbw_n_clipped_full = if (weight_type == "SBW") point_est$nclip_full else NA_integer_,
      sbw_max_abs_clipped_full = if (weight_type == "SBW") point_est$maxclip_full else NA_real_,
      sbw_n_clipped_boot_mean = if (weight_type == "SBW") mean(nclip_boot, na.rm = TRUE) else NA_real_,
      sbw_n_clipped_boot_max  = if (weight_type == "SBW") max(nclip_boot, na.rm = TRUE) else NA_integer_,
      sbw_max_abs_clipped_boot_max = if (weight_type == "SBW") max(maxclip_boot, na.rm = TRUE) else NA_real_
    ))
  }

  if (ci_method == "wald") {
    ci_log = c(loghat - z * se_boot, loghat + z * se_boot)
  } else {
    ci_log = unname(stats::quantile(boot_log[finite_rep], probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  }

  list(
    MC_fail = FALSE,
    log_est = loghat,
    est = exp(loghat),
    se_log = se_boot,
    ci_log = ci_log,
    ci = exp(ci_log),
    boot_fail_rate = fail_rate,
    boot_n_finite_reps = boot_n_finite_reps,

    # SBW clipping diagnostics
    sbw_n_clipped_full = if (weight_type == "SBW") point_est$nclip_full else NA_integer_,
    sbw_max_abs_clipped_full = if (weight_type == "SBW") point_est$maxclip_full else NA_real_,
    sbw_n_clipped_boot_mean = if (weight_type == "SBW") mean(nclip_boot, na.rm = TRUE) else NA_real_,
    sbw_n_clipped_boot_max  = if (weight_type == "SBW") max(nclip_boot, na.rm = TRUE) else NA_integer_,
    sbw_max_abs_clipped_boot_max = if (weight_type == "SBW") max(maxclip_boot, na.rm = TRUE) else NA_real_
  )
}
