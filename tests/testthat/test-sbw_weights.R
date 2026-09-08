# Regression tests for sbw_weights.R (CODE_PACKAGE_PLAN.md Step 3).
# Each test pins current behavior on a fixed toy input against hardcoded
# expected values captured from the current code, not derived independently.

test_that("get_weights_for_group_neg matches pinned closed-form weights", {
  data <- data.frame(x1 = c(1, 2, 3, 4, 5))
  X_n <- c(x1 = 2.5)

  w <- get_weights_for_group_neg(data, X_n)

  expect_equal(w, c(1.50, 1.25, 1.00, 0.75, 0.50), tolerance = 1e-8)
  # sanity: weights actually balance to the target
  expect_equal(sum(w * data$x1) / sum(w), 2.5, tolerance = 1e-8)
})

test_that(".clip_weights zeroes negatives and reports clipping diagnostics", {
  w_dirty <- c(0.5, -0.001, 0.3, -2, -1e-12, 0.9)

  out <- .clip_weights(w_dirty)

  expect_equal(out$w, c(0.5, 0, 0.3, 0, 0, 0.9), tolerance = 1e-8)
  expect_equal(out$n_clipped, 3L)
  expect_equal(out$max_abs_clipped, 2, tolerance = 1e-8)
})

test_that("get_weights_for_group_nonneg matches pinned QP weights", {
  data <- data.frame(x1 = c(1, 2, 3, 4, 5))
  X_n <- c(x1 = 2.5)

  out <- get_weights_for_group_nonneg(data, X_n)

  expect_equal(out$w, c(1.50, 1.25, 1.00, 0.75, 0.50), tolerance = 1e-8)
  expect_equal(out$n_clipped, 0L)
  expect_equal(out$max_abs_clipped, 0)
})

test_that("get_sbws_for_study accepts the closed-form solve when it's nonnegative", {
  X_subset <- data.frame(
    x1 = c(5.3, 4, 2.4, 5.3, 2.8, 7.1, 4.3, 4.2),
    x2 = c(1.6, 5.3, 4.6, 1.1, 2.1, 4.6, 5.5, 3.9)
  )
  A <- c(1, 1, 1, 1, 0, 0, 0, 0)

  out <- get_sbws_for_study(X_subset, A)

  expect_equal(
    out$w,
    c(0.986253012362326, 1.9888510002719, 0.315342655050525,
      0.709553332315246, 1.5053744504319, 1.02102218115749,
      0.465599052479489, 1.00800431593112),
    tolerance = 1e-8
  )
  expect_equal(out$n_clipped, 0L)
  expect_equal(out$max_abs_clipped, 0)
})

test_that("get_sbws_for_study falls back to nonneg QP when closed-form goes negative", {
  X_subset <- data.frame(x1 = c(1.3, 2.5, 3.2, 4.9, 1.2, 7.2, 9.0, 8.8))
  A <- c(1, 1, 1, 1, 0, 0, 0, 0)

  out <- get_sbws_for_study(X_subset, A)

  expect_equal(
    out$w,
    c(0, 0, 0.323529411764707, 3.67647058823529, 1.95368985290451,
      0.884131139366741, 0.56326352530541, 0.598915482423336),
    tolerance = 1e-8
  )
  expect_equal(out$n_clipped, 2L)
  # exact value is BLAS-dependent floating-point dust from the QP solve
  expect_lt(out$max_abs_clipped, 1e-8)
})
