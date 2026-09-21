# Shared RMST-difference estimator: E[min(T1,t*)] - E[min(T0,t*)] from
# (weighted) KM + Wald CI on the natural (not log) scale, plus its bootstrap
# wrapper. Mirrors km_ratio.R's structure exactly -- same point-estimate /
# bootstrap / fallback-on-failure shape, same SBW/IPW weight_type switch --
# but for the RMST-at-a-horizon estimand instead of the survival-ratio-at-t0
# estimand: additive (difference, not ratio), so no log-scale transform is
# needed, and the horizon is written `t_star` throughout (RMST's own horizon
# notation, kept distinct from a treatment-effect `tau` used elsewhere in
# this codebase's calling code).
#
# Added 2026-09-18 for sbw-randomization-schemes (Paper 2)'s RMST outcome
# type -- see simulations/simulation_axes.md, "Second outcome type."

#' Extract restricted mean survival time + SE at a horizon, by arm
#'
#' @param fit A `survival::survfit` object fit on `Surv(time, status) ~ A`.
#' @param t_star Horizon at which to compute the restricted mean survival time.
#' @return A list with `rmst0`, `rmst1` (RMST estimates for A=0/A=1) and
#'   `se0`, `se1` (their standard errors).
#' @keywords internal
.extract_rmst_at_tstar = function(fit, t_star) {
  tbl = summary(fit, rmean = t_star)$table
  rn = rownames(tbl)

  find_arm = function(val) {
    idx = grep(paste0("A=", val, "\\b"), rn)
    if (length(idx) == 0) idx = grep(paste0("A=", val), rn)
    if (length(idx) == 0) stop("Couldn't identify RMST strata for A=", val)
    idx[1]
  }

  i0 = find_arm(0)
  i1 = find_arm(1)

  list(
    rmst0 = tbl[i0, "rmean"], rmst1 = tbl[i1, "rmean"],
    se0 = tbl[i0, "se(rmean)"], se1 = tbl[i1, "se(rmean)"]
  )
}

#' Weighted KM restricted-mean-survival-time difference with a Wald CI
#'
#' @param time Event/censoring time.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param A Treatment indicator (0/1).
#' @param t_star Horizon at which to compute the restricted mean survival time.
#' @param alpha Significance level for the confidence interval (default 0.05).
#' @param weights Optional case weights (e.g. SBW or IPW), passed to
#'   `survival::survfit()`.
#' @return A list with `rmst0`, `rmst1`, `diff` (= RMST1 - RMST0), `se_diff`,
#'   and `ci_diff`.
#' @export
weighted_km_rmst_diff = function(time, status, A, t_star, alpha = 0.05, weights = NULL) {
  df = data.frame(time = time, status = status, A = as.integer(A))

  fit = survival::survfit(survival::Surv(time, status) ~ A, data = df, weights = weights)

  out = .extract_rmst_at_tstar(fit, t_star = t_star)

  diff = out$rmst1 - out$rmst0
  se_diff = sqrt(out$se0^2 + out$se1^2)

  z = stats::qnorm(1 - alpha / 2)
  ci_diff = c(diff - z * se_diff, diff + z * se_diff)

  list(
    rmst0 = out$rmst0, rmst1 = out$rmst1,
    diff = diff,
    se_diff = se_diff,
    ci_diff = ci_diff
  )
}

#' Bootstrap CI for a weighted-KM RMST difference (SBW or IPW)
#'
#' Bootstraps the standard error of [weighted_km_rmst_diff()]'s difference,
#' then builds either a Wald CI or a percentile CI (both on the natural,
#' additive scale -- RMST difference needs no log transform, unlike the
#' survival-ratio estimand). If the point estimate or bootstrap SE come out
#' non-finite, falls back to the unadjusted RMST difference and sets
#' `MC_fail = TRUE`; if even the unadjusted difference is undefined (a
#' double-degenerate case with no valid fallback), throws an error rather
#' than returning a nonsense finite estimate. Structurally identical to
#' [boot_km_ratio()] -- see its docs for the shared fallback logic.
#'
#' @param time Event/censoring time.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param A Treatment indicator (0/1).
#' @param X_subset Covariate data frame for the full study sample, used to
#'   fit the SBW/IPW weights.
#' @param t_star Horizon at which to compute the restricted mean survival time.
#' @param B Number of bootstrap replicates.
#' @param alpha Significance level for the confidence interval (default 0.05).
#' @param weight_type Either `"SBW"` ([get_sbws_for_study()]) or `"IPW"`
#'   ([get_ipws_for_study()]).
#' @param ci_method Either `"wald"` (default) or `"percentile"`.
#' @param ipw_use_glm Passed through to [get_ipws_for_study()] when
#'   `weight_type = "IPW"`.
#' @param seed Optional seed set at the start of the bootstrap.
#' @param verbose Currently unused; reserved for future diagnostic output.
#' @return A list with `MC_fail`, `est` (the RMST difference), `se`, `ci`,
#'   `boot_fail_rate`, `boot_n_finite_reps`, and SBW clipping diagnostics
#'   (`sbw_n_clipped_full`, `sbw_max_abs_clipped_full`,
#'   `sbw_n_clipped_boot_mean`, `sbw_n_clipped_boot_max`,
#'   `sbw_max_abs_clipped_boot_max`; all `NA` when `weight_type = "IPW"`).
#' @export
boot_km_rmst = function(time, status, A, X_subset, t_star,
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
  z = stats::qnorm(1 - alpha / 2)

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

    full = weighted_km_rmst_diff(time, status, A, t_star = t_star, alpha = alpha, weights = w_full)

    list(
      finite_rep = is.finite(full$diff),
      w_full = w_full,
      dhat = full$diff,
      nclip_full = nclip_full,
      maxclip_full = maxclip_full
    )
  }, error = function(e) {
    list(finite_rep = FALSE)
  })

  if (!isTRUE(point_est$finite_rep)) {
    # HARD FAIL, fallback to unadjusted, mark MC_fail
    unadj = tryCatch(
      weighted_km_rmst_diff(time, status, A, t_star = t_star, alpha = alpha, weights = NULL),
      error = function(e) NULL
    )
    if (is.null(unadj) || !is.finite(unadj$diff)) {
      stop("Unadjusted RMST difference failed (cannot fallback).")
    }

    return(list(
      MC_fail = TRUE,
      est = unadj$diff,
      se = NA_real_,
      ci = unadj$ci_diff,
      boot_fail_rate = NA_real_,
      boot_n_finite_reps = NA_integer_,
      sbw_n_clipped_full = NA_integer_,
      sbw_max_abs_clipped_full = NA_real_,
      sbw_n_clipped_boot_mean = NA_real_,
      sbw_n_clipped_boot_max  = NA_integer_,
      sbw_max_abs_clipped_boot_max = NA_real_
    ))
  }

  dhat = point_est$dhat

  # -------------------------
  # 2) Bootstrap attempt
  # -------------------------
  boot_diff = rep(NA_real_, B)
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

      weighted_km_rmst_diff(time_boot, status_boot, A_boot, t_star = t_star, alpha = alpha, weights = w_boot)$diff
    }, error = function(e) NA_real_)

    if (!is.finite(bootstrap_point_est)) {
      fail[b] = TRUE
      bootstrap_point_est = tryCatch(
        weighted_km_rmst_diff(time_boot, status_boot, A_boot, t_star = t_star, alpha = alpha, weights = NULL)$diff,
        error = function(e) NA_real_
      )
    }
    boot_diff[b] = bootstrap_point_est
  }

  finite_rep = is.finite(boot_diff)
  boot_n_finite_reps = sum(finite_rep)
  fail_rate = mean(fail)

  se_boot = stats::sd(boot_diff[finite_rep])

  if (!is.finite(se_boot) || se_boot <= 0) {
    unadj = tryCatch(
      weighted_km_rmst_diff(time, status, A, t_star = t_star, alpha = alpha, weights = NULL),
      error = function(e) NULL
    )
    if (is.null(unadj) || !is.finite(unadj$diff)) {
      stop("Unadjusted RMST difference failed (cannot fallback).")
    }

    return(list(
      MC_fail = TRUE,
      est = unadj$diff,
      se = NA_real_,
      ci = unadj$ci_diff,
      boot_fail_rate = fail_rate,
      boot_n_finite_reps = boot_n_finite_reps,

      sbw_n_clipped_full = if (weight_type == "SBW") point_est$nclip_full else NA_integer_,
      sbw_max_abs_clipped_full = if (weight_type == "SBW") point_est$maxclip_full else NA_real_,
      sbw_n_clipped_boot_mean = if (weight_type == "SBW") mean(nclip_boot, na.rm = TRUE) else NA_real_,
      sbw_n_clipped_boot_max  = if (weight_type == "SBW") max(nclip_boot, na.rm = TRUE) else NA_integer_,
      sbw_max_abs_clipped_boot_max = if (weight_type == "SBW") max(maxclip_boot, na.rm = TRUE) else NA_real_
    ))
  }

  ci = if (ci_method == "wald") {
    c(dhat - z * se_boot, dhat + z * se_boot)
  } else {
    unname(stats::quantile(boot_diff[finite_rep], probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
  }

  list(
    MC_fail = FALSE,
    est = dhat,
    se = se_boot,
    ci = ci,
    boot_fail_rate = fail_rate,
    boot_n_finite_reps = boot_n_finite_reps,

    sbw_n_clipped_full = if (weight_type == "SBW") point_est$nclip_full else NA_integer_,
    sbw_max_abs_clipped_full = if (weight_type == "SBW") point_est$maxclip_full else NA_real_,
    sbw_n_clipped_boot_mean = if (weight_type == "SBW") mean(nclip_boot, na.rm = TRUE) else NA_real_,
    sbw_n_clipped_boot_max  = if (weight_type == "SBW") max(nclip_boot, na.rm = TRUE) else NA_integer_,
    sbw_max_abs_clipped_boot_max = if (weight_type == "SBW") max(maxclip_boot, na.rm = TRUE) else NA_real_
  )
}
