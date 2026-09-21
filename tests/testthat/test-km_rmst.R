# Regression tests for km_rmst.R.
# Each test pins current behavior on a fixed-seed toy input against hardcoded
# expected values captured from the current code, not derived independently
# -- same convention as test-km_ratio.R.

test_that(".extract_rmst_at_tstar matches pinned survfit summary values", {
  set.seed(55)
  n <- 20
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = 0.1)
  status <- rbinom(n, 1, 0.8)
  df <- data.frame(time = time, status = status, A = A)
  fit <- survfit(Surv(time, status) ~ A, data = df)

  out <- .extract_rmst_at_tstar(fit, t_star = 5)

  expect_equal(out$rmst0, 4.49670656899439, tolerance = 1e-8)
  expect_equal(out$rmst1, 4.27074110017248, tolerance = 1e-8)
  expect_equal(out$se0, 0.474509597320943, tolerance = 1e-8)
  expect_equal(out$se1, 0.420803332313305, tolerance = 1e-8)
})

test_that("weighted_km_rmst_diff matches pinned unweighted output", {
  set.seed(123)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)

  res <- weighted_km_rmst_diff(time, status, A, t_star = 5)

  expect_equal(res$rmst0, 3.49386700756239, tolerance = 1e-8)
  expect_equal(res$rmst1, 4.68207847580722, tolerance = 1e-8)
  expect_equal(res$diff, 1.18821146824482, tolerance = 1e-8)
  expect_equal(res$se_diff, 0.456404921142446, tolerance = 1e-8)
  expect_equal(res$ci_diff, c(0.293674260438786, 2.082748676050859), tolerance = 1e-8)
})

test_that("weighted_km_rmst_diff matches pinned SBW-weighted output", {
  set.seed(123)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)
  X_subset <- data.frame(x1 = rnorm(n))
  w <- get_sbws_for_study(X_subset, A)$w

  res <- weighted_km_rmst_diff(time, status, A, t_star = 5, weights = w)

  expect_equal(res$rmst0, 3.45230749849483, tolerance = 1e-8)
  expect_equal(res$rmst1, 4.64140201458579, tolerance = 1e-8)
  expect_equal(res$diff, 1.18909451609096, tolerance = 1e-8)
  expect_equal(res$se_diff, 0.463687616333155, tolerance = 1e-8)
  expect_equal(res$ci_diff, c(0.280283488000751, 2.097905544181170), tolerance = 1e-8)
})

test_that("boot_km_rmst (SBW, Wald CI) matches pinned bootstrap output", {
  set.seed(321)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)
  X_subset <- data.frame(x1 = rnorm(n))

  res <- boot_km_rmst(time, status, A, X_subset, t_star = 5, B = 200,
                       weight_type = "SBW", ci_method = "wald", seed = 999)

  expect_false(res$MC_fail)
  expect_equal(res$est, 0.270402254828749, tolerance = 1e-8)
  expect_equal(res$se, 0.411998980478142, tolerance = 1e-8)
  expect_equal(res$ci, c(-0.537100908575631, 1.077905418233129), tolerance = 1e-8)
  expect_equal(res$boot_fail_rate, 0)
  expect_equal(res$boot_n_finite_reps, 200L)
  expect_equal(res$sbw_n_clipped_full, 0L)
  expect_equal(res$sbw_max_abs_clipped_full, 0)
  expect_equal(res$sbw_n_clipped_boot_mean, 0.015, tolerance = 1e-8)
  expect_equal(res$sbw_n_clipped_boot_max, 1L)
  # exact value is BLAS-dependent floating-point dust from the QP solve
  expect_lt(res$sbw_max_abs_clipped_boot_max, 1e-8)
})

test_that("boot_km_rmst (SBW, percentile CI) matches pinned bootstrap output", {
  set.seed(321)
  n <- 40
  A <- rep(c(0, 1), each = n / 2)
  time <- rexp(n, rate = ifelse(A == 1, 0.08, 0.12))
  status <- rbinom(n, 1, 0.85)
  X_subset <- data.frame(x1 = rnorm(n))

  res <- boot_km_rmst(time, status, A, X_subset, t_star = 5, B = 200,
                       weight_type = "SBW", ci_method = "percentile", seed = 999)

  expect_false(res$MC_fail)
  # point estimate is identical to the Wald-CI call above (same seed / same
  # bootstrap draws) -- only the CI construction differs
  expect_equal(res$est, 0.270402254828749, tolerance = 1e-8)
  expect_equal(res$ci, c(-0.413754703818217, 1.024534378633790), tolerance = 1e-8)
})

test_that("boot_km_rmst errors (rather than returning a nonsense finite estimate) on a true double-degenerate case", {
  # RMST's rmean/se(rmean) is numerically far more forgiving than the
  # survival-ratio's log-log transform (no log/division blowup at S=0 or
  # S=1), so the ratio test's "identical constant times" scenario doesn't
  # actually degenerate here. The genuine failure mode instead: a single arm
  # present at all (e.g. every subject A=1) -- survfit(Surv(...) ~ A) then
  # has only one stratum, .extract_rmst_at_tstar() can't find the other arm,
  # and neither the adjusted nor the unadjusted estimate is recoverable.
  n <- 20
  A <- rep(1, n)
  time <- rexp(n, 0.1)
  status <- rbinom(n, 1, 0.8)
  X_subset <- data.frame(x1 = rnorm(n))

  expect_error(
    boot_km_rmst(time, status, A, X_subset, t_star = 3, B = 10,
                 weight_type = "SBW", seed = 1),
    "Unadjusted RMST difference failed"
  )
})
