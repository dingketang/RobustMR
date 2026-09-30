test_that("the installed estimator agrees with the weighted moment ratio", {
  dat <- small_summary_data()
  fit <- RobustMR::mr_wald(dat)
  weights <- 1 / dat$se_gamma_ot^2
  expected <- as.numeric(crossprod(weights * dat$gamma_tr, dat$Gamma_ot) /
    crossprod(weights * dat$gamma_tr, dat$gamma_ot))

  expect_equal(fit$pe, expected, tolerance = 1e-12)
  expect_true(all(is.finite(unlist(fit))))
})

test_that("changing the external exposure scale preserves the point estimate", {
  dat <- small_summary_data()
  rescaled <- dat
  rescaled$gamma_tr <- 3 * rescaled$gamma_tr
  rescaled$se_gamma_tr <- 3 * rescaled$se_gamma_tr

  expect_equal(
    RobustMR::mr_wald(rescaled)$pe,
    RobustMR::mr_wald(dat)$pe,
    tolerance = 1e-12
  )
})

test_that("bootstrap inference is reproducible under a fixed seed", {
  dat <- small_summary_data()
  set.seed(812)
  first <- RobustMR::mr_wald_bs(dat, repit = 20L)
  set.seed(812)
  second <- RobustMR::mr_wald_bs(dat, repit = 20L)

  expect_identical(first, second)
  expect_named(first, c("pe", "lb", "ub"))
  expect_true(all(is.finite(unlist(first))))
  expect_equal(first$pe, RobustMR::mr_wald(dat)$pe)
  expect_lt(first$lb, first$pe)
  expect_gt(first$ub, first$pe)
})

test_that("the score estimator runs on a simple monotone-score example", {
  fit <- RobustMR::mr_wald_R(small_summary_data(), min_num = -1, max_num = 1)

  expect_named(fit, c("pe", "lb", "ub"))
  expect_true(all(is.finite(fit)))
  expect_true(fit[["pe"]] > 0.45 && fit[["pe"]] < 0.55)
  # This smoke test does not validate confidence-set inversion or coverage.
})
