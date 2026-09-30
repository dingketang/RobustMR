test_that("data generation has a stable schema and deterministic seed", {
  make_data <- function() {
    RobustMR::data_gen(
      seed = 42L,
      n = 150L,
      p = 7L,
      mu = 0,
      alpha_star = rep(0, 7L),
      tau0 = 0,
      gamma = seq(0.12, 0.24, length.out = 7L),
      gamma_fun = function(x) 1.3 * (x + 0.02),
      MAF = 0.3,
      beta_0 = 0.5
    )
  }
  first <- make_data()
  second <- make_data()

  expect_identical(first, second)
  expect_named(first, c("mat_all", "mat_h"))
  expect_equal(dim(first$mat_all), c(7L, 6L))
  expect_equal(dim(first$mat_h), c(7L, 4L))
  expect_named(first$mat_all, c(
    "Gamma_ot", "gamma_ot", "gamma_tr", "se_Gamma_ot", "se_gamma_tr",
    "se_gamma_ot"
  ))
  expect_named(first$mat_h, c(
    "beta.outcome", "beta.exposure", "se.outcome", "se.exposure"
  ))
  expect_true(all(is.finite(as.matrix(first$mat_all))))
  expect_equal(first$mat_h$beta.outcome, first$mat_all$Gamma_ot)
  expect_equal(first$mat_h$beta.exposure, first$mat_all$gamma_tr)
  expect_equal(first$mat_h$se.outcome, first$mat_all$se_Gamma_ot)
  expect_equal(first$mat_h$se.exposure, first$mat_all$se_gamma_tr)
})
