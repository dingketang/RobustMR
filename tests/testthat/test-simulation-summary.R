test_that("performance summaries match independently specified examples", {
  fits <- list(
    list(
      method_a = c(pe = 0.4, lb = 0.3, ub = 0.6),
      method_b = c(pe = 0.8, lb = 0.7, ub = 0.9)
    ),
    list(
      method_a = c(pe = 0.6, lb = 0.4, ub = 0.8),
      method_b = c(pe = 0.9, lb = 0.8, ub = 1.1)
    )
  )
  expected <- matrix(
    c(0, 20, 70, 100, 70, 70.7, 50, 0),
    nrow = 4L,
    dimnames = list(
      c("Bias", "RMSE", "CI length", "CI"),
      c("method_a", "method_b")
    )
  )

  expect_equal(RobustMR::process_fit_result(fits, beta_0 = 0.5), expected)
})
