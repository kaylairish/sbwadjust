# Tests for get_ipws_for_study() (R/ipw.R), the one exported function with no
# test coverage as of v0.2.0. The use_glm = FALSE path has a closed form, so
# it is checked against hand-derived values rather than pinned output; the glm
# path is cross-checked against an independently fitted stats::glm().

make_ipw_toy = function(n = 80, seed = 7) {
  set.seed(seed)
  X = data.frame(
    age = rnorm(n, 50, 10),
    bmi = rnorm(n, 27, 4)
  )
  lp = -0.5 + 0.02 * (X$age - 50) + 0.05 * (X$bmi - 27)
  A = rbinom(n, 1, plogis(lp))
  list(X = X, A = A)
}

test_that("use_glm = FALSE gives the closed-form stabilized weights", {
  d = make_ipw_toy()
  w = get_ipws_for_study(d$X, d$A, use_glm = FALSE, p_known = 0.5)

  pi_hat = mean(d$A)
  expect_equal(w[d$A == 1], rep(pi_hat / 0.5, sum(d$A == 1)))
  expect_equal(w[d$A == 0], rep((1 - pi_hat) / 0.5, sum(d$A == 0)))
  expect_length(w, length(d$A))
})

test_that("p_known and pi enter the closed form as documented", {
  d = make_ipw_toy()

  # non-0.5 known propensity, e.g. 2:1 randomization
  w = get_ipws_for_study(d$X, d$A, use_glm = FALSE, p_known = 2 / 3)
  pi_hat = mean(d$A)
  expect_equal(w[d$A == 1], rep(pi_hat / (2 / 3), sum(d$A == 1)))
  expect_equal(w[d$A == 0], rep((1 - pi_hat) / (1 - 2 / 3), sum(d$A == 0)))

  # explicit pi overrides mean(A) for the stabilization numerator only
  w_pi = get_ipws_for_study(d$X, d$A, use_glm = FALSE, p_known = 0.5, pi = 0.5)
  expect_equal(w_pi[d$A == 1], rep(0.5 / 0.5, sum(d$A == 1)))
  expect_equal(w_pi[d$A == 0], rep(0.5 / 0.5, sum(d$A == 0)))
})

test_that("the glm path matches an independently fitted propensity model", {
  d = make_ipw_toy()
  w = get_ipws_for_study(d$X, d$A)

  df = as.data.frame(d$X)
  df$A = as.integer(d$A)
  fit = stats::glm(A ~ ., family = stats::binomial(), data = df)
  ps = as.numeric(stats::predict(fit, type = "response"))
  pi_hat = mean(d$A)
  expected = ifelse(d$A == 1, pi_hat / ps, (1 - pi_hat) / (1 - ps))

  expect_equal(w, expected, tolerance = 1e-12)
  expect_true(all(is.finite(w)))
  expect_true(all(w > 0))
})
