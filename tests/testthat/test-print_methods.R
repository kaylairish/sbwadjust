# Smoke tests for the S3 print/plot methods (R/sbw_weights.R, R/sbw_estimate.R).
# These were exported but untested as of v0.2.0. They assert the identifying
# header lines and the documented invisible return, not exact formatting, so
# that cosmetic wording changes don't break the suite.

make_print_toy = function(n = 60, seed = 3) {
  set.seed(seed)
  df = data.frame(
    age = rnorm(n, 50, 10),
    bmi = rnorm(n, 27, 4),
    arm = rbinom(n, 1, 0.5)
  )
  df$Y = with(df, 0.02 * age + 0.5 * arm + rnorm(nrow(df), sd = 0.5))
  df
}

test_that("print.sbw_fit shows the header and arm counts, returning x invisibly", {
  df = make_print_toy()
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)

  expect_output(print(sbw), "<sbw_fit>", fixed = TRUE)
  expect_output(print(sbw), paste0("n:\\s+", sbw$n))
  expect_output(print(sbw), paste0(sum(sbw$treatment == 1L), " treated"))

  expect_output(out <- print(sbw))
  expect_identical(out, sbw)
})

test_that("print.sbw_fit names the treated and control levels for a character treatment", {
  df = make_print_toy()
  df$arm_chr = ifelse(df$arm == 1, "placebo", "drug")
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm_chr)

  # alphabetical: drug = 0 (control), placebo = 1 (treated)
  expect_output(print(sbw), "treated [placebo]", fixed = TRUE)
  expect_output(print(sbw), "control [drug]", fixed = TRUE)

  sbw_num = sbw_weights(~ age + bmi, data = df, treatment = arm)
  expect_false(any(grepl("[", capture.output(print(sbw_num)), fixed = TRUE)))
})

test_that("print.sbw_fit keeps a long balance formula on one line, without deparse's padding", {
  df = make_print_toy()
  for (v in c("height_cm", "weight_kg", "systolic_bp", "cholesterol")) df[[v]] = rnorm(nrow(df))
  sbw = sbw_weights(~ age + bmi + height_cm + weight_kg + systolic_bp + cholesterol,
                    data = df, treatment = arm)
  out = capture.output(print(sbw))
  expect_true(any(grepl("cholesterol", out)))
  expect_false(any(grepl("\\+\\s{2,}", out)))
})

test_that("print.summary.sbw_fit shows the balance table and ESS", {
  df = make_print_toy()
  s = summary(sbw_weights(~ age + bmi, data = df, treatment = arm))

  expect_s3_class(s, "summary.sbw_fit")
  expect_output(print(s), "Balance")
  expect_output(print(s), "Effective sample size")
  # both balance covariates are named in the printed table
  expect_output(print(s), "age")
  expect_output(print(s), "bmi")

  expect_output(out <- print(s))
  expect_identical(out, s)
})

test_that("plot.sbw_fit draws without error and returns x invisibly", {
  df = make_print_toy()
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(out <- plot(sbw))
  expect_identical(out, sbw)
})

test_that("print.sbw_estimate shows the estimand and a CI line", {
  df = make_print_toy()
  sbw = sbw_weights(~ age + bmi, data = df, treatment = arm)
  res = sbw_estimate(sbw, Y ~ 1, estimand = "ATE", B = 20, seed = 5)

  expect_output(print(res), "<sbw_estimate: ATE>", fixed = TRUE)
  expect_output(print(res), "CI:")

  expect_output(out <- print(res))
  expect_identical(out, res)
})
