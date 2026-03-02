#' @name scenario1case3
#' @title One Monte Carlo replication for Scenario 1 Case 3
#' @description
#' Runs one Monte Carlo replication under "Scenario 1, Case 3" and returns
#' point estimates and 95\% confidence intervals from multiple Mendelian
#' randomization (MR) methods. Internally, this function calls \code{\link{data_gen}},
#' then computes estimators including \code{\link{mr_wald_bs}} and \code{\link{mr_wald_qr_bs}},
#' as well as several competing MR methods.
#'
#' @param seed Integer. Seed index for this Monte Carlo replication.
#'
#' @return A named list. Each element is a numeric vector of length 3 with entries
#' \code{c(estimate, lb, ub)}:
#' \describe{
#'   \item{mr_wald}{Wald-type estimator with bootstrap-based CI (via \code{\link{mr_wald_bs}}).}
#'   \item{mr_wald_qr}{Quantile-regression Wald estimator with bootstrap-based CI (via \code{\link{mr_wald_qr_bs}}).}
#'   \item{mr_w_median}{Weighted median estimator with normal-approximation CI.}
#'   \item{mr_egger}{MR-Egger regression estimator with normal-approximation CI.}
#'   \item{mr_rap}{MR-RAPS estimator with normal-approximation CI.}
#'   \item{mr_divw}{DIVW estimator with normal-approximation CI.}
#' }
#'
#' @export
scenario1case3 <- function(seed){

  ## Sample size and number of instruments
  n <- 10000
  p <- 200

  ## Fixed RNG seed for generating gamma across Monte Carlo runs
  ## (so only the data noise varies with `seed`)
  set.seed(1234)

  ## Instrument–exposure effects (true gamma_j)
  gamma <- runif(p, min = 0.05, max = 0.1)

  ## Nonlinear transformation of gamma (used inside data generator)
  ## This induces nonlinearity / heterogeneity across instruments
  gamma_fun <- function(gamma){
    sin(gamma / 0.1 * pi / 2) * 0.2
  }

  ## No direct pleiotropic effects
  alpha_star <- rep(0, p)

  ## Generate two-sample MR data
  ## - mat_all : (Gamma_ot, gamma_ot, gamma_tr, SEs)
  ## - mat_h   : (Gamma_ot, gamma_tr, SEs)
  ##   mat_h   : (heterogenous MR data)
  data <- data_gen(seed       = seed,
                   alpha_star = alpha_star,
                   gamma_fun  = gamma_fun,
                   gamma      = gamma,
                   beta_0     = 0.5)

  ## --------------------------------------------------------------------------
  ## Proposed methods
  ## --------------------------------------------------------------------------

  ## Wald-type estimator with bootstrap CI
  fit1 <- mr_wald_bs(data$mat_all)

  ## Quantile-regression-based Wald estimator with bootstrap CI
  fit2 <- mr_wald_qr_bs(data$mat_all)

  ## Outcome–exposure summary data for standard MR methods
  mat_h <- data$mat_h

  ## --------------------------------------------------------------------------
  ## Competing MR estimators
  ## --------------------------------------------------------------------------

  ## Weighted median estimator
  fit3 <- mr_weighted_median(
    b_exp  = mat_h$beta.exposure,
    b_out  = mat_h$beta.outcome,
    se_exp = mat_h$se.exposure,
    se_out = mat_h$se.outcome
  )

  ## MR-RAPS with Tukey loss (robust to weak instruments / outliers)
  fit4 <- mr.raps(mat_h,
                  loss.function  = "tukey",
                  diagnostics    = FALSE,
                  over.dispersion = TRUE)

  ## MR-Egger regression
  fit5 <- mr_egger_regression(
    b_exp  = mat_h$beta.exposure,
    b_out  = mat_h$beta.outcome,
    se_exp = mat_h$se.exposure,
    se_out = mat_h$se.outcome
  )

  ## DIVW estimator
  fit6 <- mr.divw::mr.divw(
    beta.exposure = mat_h$beta.exposure,
    beta.outcome  = mat_h$beta.outcome,
    se.exposure   = mat_h$se.exposure,
    se.outcome    = mat_h$se.outcome
  )

  ## --------------------------------------------------------------------------
  ## Collect point estimates and 95% confidence intervals
  ## Format: c(estimate, lower, upper)
  ## --------------------------------------------------------------------------

  mr_wald     <- c(fit1$pe,
                   lb = fit1$lb,
                   ub = fit1$ub)

  mr_wald_qr  <- c(fit2$pe,
                   lb = fit2$lb,
                   ub = fit2$ub)

  mr_w_median <- c(fit3$b,
                   lb = fit3$b - fit3$se * 1.96,
                   ub = fit3$b + fit3$se * 1.96)

  mr_rap      <- c(fit4$beta.hat,
                   fit4$beta.hat - fit4$beta.se * 1.96,
                   fit4$beta.hat + fit4$beta.se * 1.96)

  mr_egger    <- c(fit5$b,
                   fit5$b - fit5$se * 1.96,
                   fit5$b + fit5$se * 1.96)

  mr_divw     <- c(fit6$beta.hat,
                   fit6$beta.hat - fit6$beta.se * 1.96,
                   fit6$beta.hat + fit6$beta.se * 1.96)

  ## Return all estimators in a named list
  return(list(
    mr_wald     = mr_wald,
    mr_wald_qr  = mr_wald_qr,
    mr_w_median = mr_w_median,
    mr_egger    = mr_egger,
    mr_rap      = mr_rap,
    mr_divw     = mr_divw
  ))
}


#' @name process_fit_result
#' @title Summarize Monte Carlo performance metrics
#' @description
#' Post-processes a list of Monte Carlo outputs (as returned by
#' \code{\link{scenario1case3}}) to compute performance metrics for each estimator,
#' including relative bias, relative RMSE, relative CI length, and empirical
#' coverage probability.
#'
#' @param fitlist A list of Monte Carlo replications. Each element should be the
#'   output of \code{\link{scenario1case3}}, i.e., a named list of length-3 numeric
#'   vectors \code{c(estimate, lb, ub)}.
#' @param beta_0 Numeric. True causal effect used to compute relative metrics.
#'
#' @return A numeric matrix with 4 rows and one column per estimator. Rows are:
#' \describe{
#'   \item{Bias}{Relative bias in percent.}
#'   \item{RMSE}{Relative RMSE in percent.}
#'   \item{CI length}{Relative mean CI length in percent.}
#'   \item{CI}{Empirical coverage in percent.}
#' }
#' Column names correspond to estimator names.
#'
#' @export
process_fit_result <- function(fitlist, beta_0){

  ## Number of estimators
  n_estimator <- length(fitlist[[1]])

  ## Estimator names (consistent across Monte Carlo runs)
  namelist <- names(fitlist[[1]])

  estimators <- list()

  ## Performance metrics
  CI = CIlength = RMSE = Bias = rep(0, n_estimator)

  for (i in 1:n_estimator) {

    ## Stack Monte Carlo results:
    ## rows = replications, columns = (estimate, lb, ub)
    estimators[[i]] <-
      as.matrix(do.call(rbind, lapply(fitlist, `[[`, namelist[[i]])))

    ## Relative bias (%)
    Bias[i] <-
      mean(estimators[[i]][, 1] - beta_0) / beta_0 * 100

    ## Relative RMSE (%)
    RMSE[i] <-
      sqrt(mean((estimators[[i]][, 1] - beta_0)^2)) / beta_0 * 100

    ## Relative CI length (%)
    CIlength[i] <-
      mean(estimators[[i]][, 3] - estimators[[i]][, 2]) / beta_0 * 100

    ## Empirical coverage probability
    CI[i] <-
      mean((estimators[[i]][, 3] > beta_0) &
             (estimators[[i]][, 2] < beta_0))
  }

  ## Assemble summary table
  dat <- matrix(0, nrow = 4, ncol = n_estimator)
  dat[1, ] <- round(Bias, 1)
  dat[2, ] <- round(RMSE, 1)
  dat[3, ] <- round(CIlength, 1)
  dat[4, ] <- round(CI * 100, 1)

  rownames(dat) <- c("Bias", "RMSE", "CI length", "CI")
  colnames(dat) <- namelist

  dat
}

