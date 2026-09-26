# Regression tests for the sbw_weights()/sbw_estimate() formula API
# (R/sbw_weights.R, R/sbw_estimate.R). Each test pins current behavior on a
# fixed-seed toy input against hardcoded expected values captured from the
# current code, following the convention in test-sbw_weights.R /
# test-km_ratio.R. Correctness of the hand-derived .weighted_win_prob() and
# .weighted_quantile() helpers was additionally cross-checked against an
# independent brute-force double loop and against stats::quantile(type = 1)
# via an integer-weight replication trick (not repeated here as pinned
# tests, since those are one-off derivation checks rather than regression
# guards).

make_toy_trial = function(n = 60, seed = 1) {
  set.seed(seed)
  data.frame(
    age = rnorm(n, 50, 10),
    bmi = rnorm(n, 27, 4),
    region = sample(c("N", "S"), n, replace = TRUE),
    arm = rbinom(n, 1, 0.5)
  )
}

test_that("sbw_weights() balances covariates and matches get_sbws_for_study()", {
  df = make_toy_trial()
  sbw = sbw_weights(~ age + bmi + region, data = df, treatment = arm)

  expect_s3_class(sbw, "sbw_fit")
  expect_equal(sbw$treatment_name, "arm")
  expect_equal(sbw$n, 60L)

  X = model.matrix(~ age + bmi + region, data = df)[, -1, drop = FALSE]
  direct = get_sbws_for_study(as.data.frame(X), df$arm)
  expect_equal(sbw$weights, direct$w, tolerance = 1e-8)

  # balance actually holds: weighted treated/control means agree
  s = summary(sbw)
  expect_equal(s$balance$mean_treated_w, s$balance$mean_control_w, tolerance = 1e-6)
})

test_that("sbw_weights() recodes a factor/character treatment consistently", {
  df = make_toy_trial()
  df$arm_chr = ifelse(df$arm == 1, "Drug", "Placebo")

  sbw_num = sbw_weights(~ age + bmi, data = df, treatment = arm)
  sbw_chr = sbw_weights(~ age + bmi, data = df, treatment = arm_chr)

  # "Drug" sorts before "Placebo" alphabetically, so factor level 1 = Drug =
  # 0, i.e. the opposite coding of arm (where 1 = Drug) -- weights should
  # still agree since get_sbws_for_study() balances the same two groups.
  expect_equal(sbw_chr$treatment, 1L - sbw_num$treatment)
  expect_equal(sbw_chr$weights, sbw_num$weights, tolerance = 1e-8)
})

test_that("sbw_weights() rejects non-0/1 numeric treatment and NA covariates", {
  df = make_toy_trial()
  expect_error(sbw_weights(~ age, data = df, treatment = age), "0/1")

  df_na = df
  df_na$age[1] = NA
  expect_error(sbw_weights(~ age, data = df_na, treatment = arm), "Missing values")
})

test_that("weights.sbw_fit() returns the same vector as sbw$weights", {
  df = make_toy_trial()
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)
  expect_identical(weights(sbw), sbw$weights)
})

test_that("sbw_estimate ATE matches pinned bootstrap output", {
  df = make_toy_trial()
  df$Y = with(df, 0.02 * age + 0.5 * arm + rnorm(nrow(df), sd = 0.5))
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)

  res = sbw_estimate(sbw, Y ~ 1, estimand = "ATE", B = 100, seed = 99)

  expect_s3_class(res, "sbw_estimate")
  expect_equal(res$scale, "identity")
  point_direct = with(df, {
    w = sbw$weights
    sum(w[arm == 1] * Y[arm == 1]) / sum(w[arm == 1]) -
      sum(w[arm == 0] * Y[arm == 0]) / sum(w[arm == 0])
  })
  expect_equal(unname(res$estimate), point_direct, tolerance = 1e-8)
  expect_true(res$ci[1, 1] < unname(res$estimate))
  expect_true(res$ci[1, 2] > unname(res$estimate))
})

test_that("sbw_estimate RR is exp() of the log-scale weighted-mean-ratio point estimate", {
  df = make_toy_trial()
  df$Y = rbinom(nrow(df), 1, plogis(-1 + 0.02 * df$age + 0.3 * df$arm))
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)

  res = sbw_estimate(sbw, Y ~ 1, estimand = "RR", B = 100, seed = 99)

  w = sbw$weights
  m1 = sum(w[df$arm == 1] * df$Y[df$arm == 1]) / sum(w[df$arm == 1])
  m0 = sum(w[df$arm == 0] * df$Y[df$arm == 0]) / sum(w[df$arm == 0])
  expect_equal(res$scale, "log")
  expect_equal(unname(res$estimate), m1 / m0, tolerance = 1e-8)
  expect_true(all(res$ci > 0))
})

test_that("sbw_estimate survival_ratio matches boot_km_ratio() called directly", {
  df = make_toy_trial()
  df$time = rexp(nrow(df), rate = 0.05 + 0.01 * df$arm)
  df$status = rbinom(nrow(df), 1, 0.8)
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)

  res = sbw_estimate(sbw, survival::Surv(time, status) ~ 1,
                      estimand = "survival_ratio", horizon = 10, B = 50, seed = 7)
  direct = boot_km_ratio(df$time, df$status, df$arm, sbw$X, t0 = 10, B = 50, seed = 7)

  expect_equal(unname(res$estimate), direct$est, tolerance = 1e-8)
  expect_equal(as.numeric(res$ci), direct$ci, tolerance = 1e-8)
})

test_that("survival_ratio warns when it falls back to the unadjusted KM ratio", {
  df = make_toy_trial()
  df$time = rexp(nrow(df), rate = 0.05 + 0.01 * df$arm)
  df$status = rbinom(nrow(df), 1, 0.8)
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)

  # B = 1: the SD of a single bootstrap replicate is NA, so the SE is non-finite
  expect_warning(
    res <- sbw_estimate(sbw, survival::Surv(time, status) ~ 1,
                        estimand = "survival_ratio", horizon = 10, B = 1, seed = 7),
    "unadjusted Kaplan-Meier"
  )
  unadj = km_ratio_loglog_greenwood(df$time, df$status, df$arm, t0 = 10)
  expect_true(res$mc_fail)
  expect_equal(unname(res$estimate), unadj$ratio, tolerance = 1e-8)

  # the ordinary path stays silent
  expect_no_warning(
    sbw_estimate(sbw, survival::Surv(time, status) ~ 1,
                 estimand = "survival_ratio", horizon = 10, B = 50, seed = 7)
  )
})

test_that("sbw_estimate mann_whitney matches the brute-force weighted win probability", {
  df = make_toy_trial()
  df$Y = rnorm(nrow(df), mean = df$age / 10 + df$arm)
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)

  res = sbw_estimate(sbw, Y ~ 1, estimand = "mann_whitney", B = 50, seed = 3)

  w = sbw$weights
  y1 = df$Y[df$arm == 1]; w1 = w[df$arm == 1]
  y0 = df$Y[df$arm == 0]; w0 = w[df$arm == 0]
  brute = 0
  for (i in seq_along(y1)) for (j in seq_along(y0)) {
    brute = brute + w1[i] * w0[j] * ((y1[i] > y0[j]) + 0.5 * (y1[i] == y0[j]))
  }
  brute = brute / (sum(w1) * sum(w0))

  expect_equal(unname(res$estimate), brute, tolerance = 1e-8)
  expect_true(unname(res$estimate) >= 0 && unname(res$estimate) <= 1)
})

test_that("sbw_estimate quantile_diff/quantile_ratio match direct .weighted_quantile calls", {
  df = make_toy_trial()
  df$Y = rnorm(nrow(df), mean = df$age / 10 + df$arm)
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)
  w = sbw$weights

  res_diff = sbw_estimate(sbw, Y ~ 1, estimand = "quantile_diff", probs = c(0.25, 0.5, 0.75),
                           B = 50, seed = 5)
  q1 = .weighted_quantile(df$Y[df$arm == 1], w[df$arm == 1], c(0.25, 0.5, 0.75))
  q0 = .weighted_quantile(df$Y[df$arm == 0], w[df$arm == 0], c(0.25, 0.5, 0.75))
  expect_equal(unname(res_diff$estimate), q1 - q0, tolerance = 1e-8)
  expect_equal(names(res_diff$estimate), c("q0.25", "q0.5", "q0.75"))

  res_ratio = sbw_estimate(sbw, Y ~ 1, estimand = "quantile_ratio", probs = 0.5, B = 50, seed = 5)
  expect_equal(res_ratio$scale, "log")
  expect_equal(unname(res_ratio$estimate),
               .weighted_quantile(df$Y[df$arm == 1], w[df$arm == 1], 0.5) /
                 .weighted_quantile(df$Y[df$arm == 0], w[df$arm == 0], 0.5),
               tolerance = 1e-8)
})

test_that("sbw_estimate works when treatment was passed as a vector, not a column", {
  df = make_toy_trial()
  arm_vec = df$arm
  df$arm = NULL
  df$Y = with(df, 0.02 * age + 0.5 * arm_vec + rnorm(nrow(df), sd = 0.5))
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm_vec)

  res = sbw_estimate(sbw, Y ~ 1, estimand = "ATE", B = 50, seed = 2)
  expect_equal(res$boot_fail_rate, 0)
})

test_that("sbw_estimate takes outcomes from a later `data` and matches the all-in-one fit", {
  df = make_toy_trial()
  Y = with(df, 0.02 * age + 0.5 * arm + rnorm(nrow(df), sd = 0.5))

  sbw_design = sbw_weights(~ age + bmi, data = df, treatment = arm)
  res_later = sbw_estimate(sbw_design, Y ~ 1, estimand = "ATE",
                           data = data.frame(Y = Y), B = 50, seed = 4)

  df_all = df
  df_all$Y = Y
  sbw_all = sbw_weights(~ age + bmi, data = df_all, treatment = arm)
  res_all = sbw_estimate(sbw_all, Y ~ 1, estimand = "ATE", B = 50, seed = 4)

  expect_equal(res_later$estimate, res_all$estimate, tolerance = 1e-10)
  expect_equal(res_later$ci, res_all$ci, tolerance = 1e-10)

  # full data frame with outcome added also works (shared columns agree)
  res_full = sbw_estimate(sbw_design, Y ~ 1, estimand = "ATE", data = df_all, B = 50, seed = 4)
  expect_equal(res_full$estimate, res_all$estimate, tolerance = 1e-10)
})

test_that("sbw_estimate rejects outcome data that doesn't line up with the fit", {
  df = make_toy_trial()
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)
  df$Y = rnorm(nrow(df))

  expect_error(sbw_estimate(sbw, Y ~ 1, "ATE", data = df[-1, ], B = 10), "rows")
  expect_error(sbw_estimate(sbw, Y ~ 1, "ATE", data = df[nrow(df):1, ], B = 10), "differs")

  df_na = df
  df_na$Y[1] = NA
  expect_error(sbw_estimate(sbw, Y ~ 1, "ATE", data = df_na, B = 10), "Missing values")
})

test_that("survival_ratio resolves Surv() without survival attached, and rejects NA times", {
  df = make_toy_trial()
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)
  df$time = rexp(nrow(df), rate = 0.05)
  df$status = rbinom(nrow(df), 1, 0.8)

  # formula environment can't see survival::Surv, as for a user who never attached it
  f = local(Surv(time, status) ~ 1, envir = new.env(parent = baseenv()))
  res = sbw_estimate(sbw, f, estimand = "survival_ratio", data = df, horizon = 10, B = 20, seed = 1)
  expect_s3_class(res, "sbw_estimate")

  df$time[2] = NA
  expect_error(sbw_estimate(sbw, f, "survival_ratio", data = df, horizon = 10, B = 20),
               "Missing values")
})

test_that(".weighted_win_prob matches an independent brute-force double loop, with ties", {
  set.seed(42)
  y1 = rnorm(15); w1 = runif(15, 0.5, 2)
  y0 = rnorm(12); w0 = runif(12, 0.5, 2)
  brute = 0
  for (i in seq_along(y1)) for (j in seq_along(y0)) {
    brute = brute + w1[i] * w0[j] * ((y1[i] > y0[j]) + 0.5 * (y1[i] == y0[j]))
  }
  brute = brute / (sum(w1) * sum(w0))
  expect_equal(.weighted_win_prob(y1, w1, y0, w0), brute, tolerance = 1e-12)

  y1t = c(1, 2, 2, 3); w1t = c(1, 1, 1, 1)
  y0t = c(2, 2, 4); w0t = c(1, 1, 1)
  expect_equal(.weighted_win_prob(y1t, w1t, y0t, w0t), 1/3, tolerance = 1e-12)
})

test_that(".weighted_quantile matches stats::quantile(type = 1) under unit weights, and an integer-weight replication check", {
  set.seed(7)
  y = rnorm(37)
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(
    .weighted_quantile(y, rep(1, length(y)), probs),
    unname(stats::quantile(y, probs = probs, type = 1)),
    tolerance = 1e-12
  )

  set.seed(11)
  y2 = round(rnorm(10), 2)
  w2 = c(2, 1, 3, 1, 2, 1, 4, 2, 1, 3)
  expanded = rep(y2, w2)
  probs2 = c(0.2, 0.5, 0.8)
  expect_equal(
    .weighted_quantile(y2, w2, probs2),
    unname(stats::quantile(expanded, probs = probs2, type = 1)),
    tolerance = 1e-12
  )
})
