# Regression tests for km_ratio.R (CODE_PACKAGE_PLAN.md Step 3).
# Each test pins current behavior on a fixed-seed toy input against hardcoded
# expected values captured from the current code, not derived independently.

test_that(".extract_surv_at_t0 matches pinned survfit summary values", {
  set.seed(55)
  n <- 20
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = 0.1)
  status <- rbinom(n, 1, 0.8)
  df <- data.frame(time = time, status = status, A = A)
  fit <- survfit(Surv(time, status) ~ A, data = df)

  out <- .extract_surv_at_t0(fit, t0 = 5)

  expect_equal(out$S0, 0.888888888888889, tolerance = 1e-8)
  expect_equal(out$S1, 0.685714285714286, tolerance = 1e-8)
  expect_equal(out$se0, 0.104756560175785, tolerance = 1e-8)
  expect_equal(out$se1, 0.151494017432321, tolerance = 1e-8)
})

test_that("km_ratio_loglog_greenwood matches pinned unweighted output", {
  set.seed(123)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)

  res <- km_ratio_loglog_greenwood(time, status, A, t0 = 5)

  expect_equal(res$S0, 0.373333333333333, tolerance = 1e-8)
  expect_equal(res$S1, 0.9, tolerance = 1e-8)
  expect_equal(res$ratio, 2.41071428571429, tolerance = 1e-8)
  expect_equal(res$log_ratio, 0.87992308770328, tolerance = 1e-8)
  expect_equal(res$se_log, 0.306995165672086, tolerance = 1e-8)
  expect_equal(res$ci_ratio, c(1.32078151676736, 4.40007926638074), tolerance = 1e-8)
  expect_equal(res$ci_log, c(0.278223619558084, 1.48162255584848), tolerance = 1e-8)
})

test_that("km_ratio_loglog_greenwood matches pinned SBW-weighted output", {
  set.seed(123)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)
  X_subset <- data.frame(x1 = rnorm(n))
  w <- get_sbws_for_study(X_subset, A)$w

  res <- km_ratio_loglog_greenwood(time, status, A, t0 = 5, weights = w)

  expect_equal(res$S0, 0.35943176127448, tolerance = 1e-8)
  expect_equal(res$S1, 0.886500789461794, tolerance = 1e-8)
  expect_equal(res$ratio, 2.46639525210132, tolerance = 1e-8)
  expect_equal(res$log_ratio, 0.90275767255605, tolerance = 1e-8)
  expect_equal(res$se_log, 0.319306898113021, tolerance = 1e-8)
  expect_equal(res$ci_ratio, c(1.31907093575437, 4.6116591418255), tolerance = 1e-8)
  expect_equal(res$ci_log, c(0.276927652239329, 1.52858769287277), tolerance = 1e-8)
})

test_that("boot_km_ratio (SBW, Wald CI) matches pinned bootstrap output", {
  set.seed(321)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)
  X_subset <- data.frame(x1 = rnorm(n))

  res <- boot_km_ratio(time, status, A, X_subset, t0 = 5, B = 200,
                        weight_type = "SBW", ci_method = "wald", seed = 999)

  expect_false(res$MC_fail)
  expect_equal(res$log_est, 0.00371710267214298, tolerance = 1e-8)
  expect_equal(res$est, 1.00372401966602, tolerance = 1e-8)
  expect_equal(res$se_log, 0.276168969531657, tolerance = 1e-8)
  expect_equal(res$ci_log, c(-0.537564131257444, 0.54499833660173), tolerance = 1e-8)
  expect_equal(res$ci, c(0.584169480887049, 1.72460551366822), tolerance = 1e-8)
  expect_equal(res$boot_fail_rate, 0)
  expect_equal(res$boot_n_finite_reps, 200L)
  expect_equal(res$sbw_n_clipped_full, 0L)
  expect_equal(res$sbw_max_abs_clipped_full, 0)
  expect_equal(res$sbw_n_clipped_boot_mean, 0.015, tolerance = 1e-8)
  expect_equal(res$sbw_n_clipped_boot_max, 1L)
  # exact value is BLAS-dependent floating-point dust from the QP solve
  expect_lt(res$sbw_max_abs_clipped_boot_max, 1e-8)
})

test_that("boot_km_ratio (SBW, percentile CI) matches pinned bootstrap output", {
  set.seed(321)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)
  X_subset <- data.frame(x1 = rnorm(n))

  res <- boot_km_ratio(time, status, A, X_subset, t0 = 5, B = 200,
                        weight_type = "SBW", ci_method = "percentile", seed = 999)

  expect_false(res$MC_fail)
  # point estimate is identical to the Wald-CI call above (same seed / same
  # bootstrap draws) -- only the CI construction differs
  expect_equal(res$log_est, 0.00371710267214298, tolerance = 1e-8)
  expect_equal(res$est, 1.00372401966602, tolerance = 1e-8)
  expect_equal(res$ci_log, c(-0.454077087202489, 0.557809479152537), tolerance = 1e-8)
  expect_equal(res$ci, c(0.635033778381594, 1.74684181281927), tolerance = 1e-8)
})

test_that("boot_km_ratio errors (rather than returning a nonsense finite estimate) on a true double-degenerate case", {
  n <- 20
  A <- rep(c(0, 1), each = n / 2)
  time <- c(rep(1, n / 2), rep(2, n / 2))
  status <- rep(1, n)
  X_subset <- data.frame(x1 = rnorm(n))

  expect_error(
    boot_km_ratio(time, status, A, X_subset, t0 = 3, B = 50,
                  weight_type = "SBW", seed = 1),
    "Unadjusted KM ratio failed"
  )
})
