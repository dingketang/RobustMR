#' @import quantreg
#' @import MASS
#' @import TwoSampleMR
#' @import mr.raps

#' @name mr_wald
#' @title MR Wald-type estimator via two weighted no-intercept regressions
#' @description
#' Computes a Wald-type Mendelian randomization (MR) estimate using two
#' weighted linear regressions through the origin. The point estimate is the
#' ratio of the two fitted slopes.
#'
#' @param data_mat A data.frame containing columns Gamma_ot, gamma_ot, gamma_tr,
#'   se_Gamma_ot, and se_gamma_ot.
#'
#' @return A list with elements pe, lb, and ub.
#' @export
mr_wald <- function(data_mat){
  ## Fit 2: weighted regression for outcome association on instrument strength
  ## weights = 1 / Var(Gamma_ot)
  fit2_weights <- 1 / data_mat$se_Gamma_ot^2
  fit2 <- summary(lm(data_mat$Gamma_ot ~ data_mat$gamma_tr + 0,
                     weights = fit2_weights))

  ## Fit 1: weighted regression for (other) exposure association on instrument strength
  ## weights = 1 / Var(gamma_ot)
  fit1_weights <- 1 / data_mat$se_gamma_ot^2
  fit1 <- summary(lm(data_mat$gamma_ot ~ data_mat$gamma_tr + 0,
                     weights = fit1_weights))

  ## Ratio-of-slopes Wald estimate
  pe <- fit2$coefficients[1] / fit1$coefficients[1]

  ## Your current SE: uses SE of fit2 slope divided by fit1 slope
  ## (Note: this ignores uncertainty in fit1 slope.)
  se <- fit2$coefficients[2] / fit1$coefficients[1]

  return(list(
    pe   = pe,
    lb   = pe - se * 1.96,
    ub   = pe + se * 1.96,
    bhat = fit2$coefficients[1]  # slope from outcome regression
  ))
}


#' @name mr_wald_bs
#' @title Bootstrap CI for MR Wald estimator
#' @description
#' Computes a bootstrap-based normal-approximation confidence interval for
#' \code{\link{mr_wald}} by resampling SNP rows of \code{data_mat} with replacement.
#'
#' @param data_mat A data.frame containing the columns required by \code{\link{mr_wald}}.
#' @param repit Integer. Number of bootstrap replicates. Default is 500.
#'
#' @return A list with elements \code{pe}, \code{lb}, and \code{ub}.
#'
#' @export
mr_wald_bs <- function(data_mat, repit = 500){

  result <- rep(0, repit)   # store bootstrap estimates
  p <- nrow(data_mat)       # number of SNPs (rows)

  for (i in 1:repit) {
    ## resample SNPs with replacement
    index <- sample(1:p, replace = TRUE)
    result[i] <- mr_wald(data_mat[index, ])$pe
  }

  pe <- mr_wald(data_mat)$pe
  sdCI   = sd(result)
  return(list(pe= pe,lb =pe -sdCI*1.96 ,ub = pe + sdCI*1.96 ))

}


#' @name mr_wald_qr
#' @title Quantile-regression analogue of MR Wald estimator (median regression)
#' @description
#' Computes a Wald-type Mendelian randomization (MR) estimate using median
#' regression (quantile regression with \eqn{\tau = 0.5}) through the origin for
#' both outcome and exposure associations, and takes the ratio of fitted slopes.
#'
#' @param data_mat A data.frame containing columns \code{Gamma_ot}, \code{gamma_ot},
#'   \code{gamma_tr}, \code{se_Gamma_ot}, and \code{se_gamma_ot}.
#'
#' @return A list with element:
#' \describe{
#'   \item{pe}{Point estimate of the causal effect (ratio of median-regression slopes).}
#' }
#'
#' @export
mr_wald_qr <- function(data_mat){

  ## Median regression for outcome association
  fit2_weights <- 1 / data_mat$se_Gamma_ot
  fit2 <- rq(data_mat$Gamma_ot ~ data_mat$gamma_tr + 0,
             weights = fit2_weights, tau = 0.5)

  ## Median regression for gamma_ot on gamma_tr
  fit1_weights <- 1 / data_mat$se_gamma_ot
  fit1 <- rq(data_mat$gamma_ot ~ data_mat$gamma_tr + 0,
             weights = fit1_weights, tau = 0.5)

  ## Ratio of slopes
  pe <- coef(fit2, tau = 0.5)[1] / coef(fit1, tau = 0.5)[1]

  return(list(pe = pe))
}


#' @name mr_wald_qr_bs
#' @title Bootstrap CI for quantile-regression Wald estimator
#' @description
#' Computes a bootstrap-based normal-approximation confidence interval for
#' \code{\link{mr_wald_qr}} by resampling SNP rows of \code{data_mat} with replacement.
#'
#' @param data_mat A data.frame containing the columns required by \code{\link{mr_wald_qr}}.
#' @param repit Integer. Number of bootstrap replicates. Default is 500.
#'
#' @return A list with elements \code{pe}, \code{lb}, and \code{ub}.
#'
#' @export
mr_wald_qr_bs <- function(data_mat, repit = 500){

  result <- rep(0, repit)
  p <- nrow(data_mat)

  for (i in 1:repit) {
    index <- sample(1:p, replace = TRUE)
    result[i] <- mr_wald_qr(data_mat[index, ])$pe
  }

  pe <- mr_wald_qr(data_mat)$pe
  sdCI   = sd(result)
  return(list(pe= pe,lb =pe -sdCI*1.96 ,ub = pe + sdCI*1.96 )
  )
}


#' @name data_gen
#' @title Generate two-sample Mendelian randomization summary statistics
#' @description
#' Generates simulated genotype data and traits in two independent samples
#' (an "outcome sample" and a "treatment sample"), then computes SNP-wise marginal
#' associations and standard errors via simple linear regressions. Returns two
#' data.frames: one in a TwoSampleMR-style format and one containing additional
#' summary statistics used by methods in this package.
#'
#' @param seed Integer or \code{NULL}. Random seed. If \code{NULL}, a seed is generated internally.
#' @param n Integer. Sample size (used for both the outcome and treatment samples). Default is \code{10000}.
#' @param p Integer. Number of SNPs. Default is \code{200}.
#' @param mu Numeric. Mean used when generating the direct-effect vector \code{alpha}. Default is \code{0}.
#' @param alpha_star Numeric vector or \code{NULL}. Currently included as an input but not used by the function.
#' @param tau0 Numeric. Standard deviation used when generating \code{alpha}. Default is \code{0}.
#' @param gamma Numeric vector of length \code{p}. SNP-exposure effects in the treatment sample.
#' @param gamma_fun Function applied to \code{gamma} when generating \code{D} in the outcome sample.
#' @param MAF Numeric in (0, 1). Minor allele frequency used for genotype simulation. Default is \code{0.3}.
#' @param beta_0 Numeric. True causal effect used in the outcome model.
#'
#' @return A list with elements:
#' \describe{
#'   \item{mat_all}{A data.frame with columns \code{Gamma_ot}, \code{gamma_ot}, \code{gamma_tr},
#'     \code{se_Gamma_ot}, \code{se_gamma_tr}, and \code{se_gamma_ot}.}
#'   \item{mat_h}{A data.frame with columns \code{beta.outcome}, \code{beta.exposure},
#'     \code{se.outcome}, and \code{se.exposure}.}
#' }
#'
#' @export
data_gen <- function(seed       = NULL,
                     n          = 10000,
                     p          = 200,
                     mu         = 0,
                     alpha_star = NULL,
                     tau0 = 0,
                     gamma,
                     gamma_fun,
                     MAF        = 0.3,
                     beta_0){

  ## Seed handling (deterministic but different across runs)
  if (is.null(seed))
    seed <- floor(runif(1) * 10000)
  set.seed(seed * 2025)

  ## Outcome sample genotype matrix: n x p, genotype in {0,1,2}
  Z <- matrix(rbinom(n * p, 2, MAF), ncol = p)

  ## Unobserved confounder
  U <- rnorm(n)

  ## Exposure in outcome sample
  ## NOTE: gamma_fun(gamma) must exist; otherwise replace with gamma
  D <- Z %*% gamma_fun(gamma) + U + rnorm(n)

  ## Outcome in outcome sample
  ##
  alpha = rnorm(p,mean  = mu , sd = tau0)

  Y <- beta_0 * D + U + rnorm(n) + Z %*% alpha

  ## Treatment sample genotypes and exposure
  Z_new <- matrix(rbinom(n * p, 2, MAF), ncol = p)
  U <- rnorm(n)
  D_new <- Z_new %*% gamma + U + rnorm(n)

  ## Helper: SNP -> (estimate, SE) for association with D (outcome sample)
  get_gamma_ot <- function(V){
    summary(lm(D ~ V))$coefficients[2, 1:2]
  }
  gamma_ot <- apply(Z, 2, get_gamma_ot)
  se_gamma_ot <- gamma_ot[2, ]
  gamma_ot <- gamma_ot[1, ]

  ## Helper: SNP -> (estimate, SE) for association with Y (outcome sample)
  get_Gamma_ot <- function(V){
    summary(lm(Y ~ V))$coefficients[2, 1:2]
  }
  Gamma_ot <- apply(Z, 2, get_Gamma_ot)
  se_Gamma_ot <- Gamma_ot[2, ]
  Gamma_ot <- Gamma_ot[1, ]

  ## Helper: SNP -> (estimate, SE) for association with D_new (treatment sample)
  get_gamma_tr <- function(V){
    summary(lm(D_new ~ V))$coefficients[2, 1:2]
  }
  gamma_tr <- apply(Z_new, 2, get_gamma_tr)
  se_gamma_tr <- gamma_tr[2, ]
  gamma_tr <- gamma_tr[1, ]

  ## A TwoSampleMR-style data frame for outcome vs exposure summary stats
  mat_h <- data.frame(beta.outcome  = Gamma_ot,
                      beta.exposure = gamma_tr,
                      se.outcome    = se_Gamma_ot,
                      se.exposure   = se_gamma_tr)



  ## Another data.frame that report additional summary statistics
  mat_all = data.frame(Gamma_ot    = Gamma_ot,
                       gamma_ot    = gamma_ot,
                       gamma_tr    = gamma_tr,
                       se_Gamma_ot = se_Gamma_ot,
                       se_gamma_tr = se_gamma_tr,
                       se_gamma_ot = se_gamma_ot)
  return(list(
    mat_all = mat_all,
    mat_h   = mat_h
  ))
}


