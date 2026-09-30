
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
#' @details This preserves the source slope-based interval, which ignores
#' uncertainty in the denominator slope and covariance between slopes. A
#' zero or negative denominator can produce non-finite or reversed bounds.
#' Use \code{\link{mr_wald_bs}} for the separately implemented SNP-pairs
#' bootstrap interval.
#'
#' @return A list with elements \code{pe}, \code{lb}, \code{ub}, and
#'   \code{bhat} (the outcome-regression slope).
#'
#' @examples
#' x <- seq(0.08, 0.24, length.out = 24)
#' dat <- data.frame(gamma_tr = x, gamma_ot = 1.3 * x,
#'   Gamma_ot = 0.5 * 1.3 * x + 0.008 * sin(seq_along(x)),
#'   se_gamma_ot = 0.03, se_gamma_tr = 0.02, se_Gamma_ot = 0.04)
#' mr_wald(dat)
#' @export
mr_wald <- function(data_mat){
  ## Fit 2: weighted regression for outcome association on instrument strength
  ## Both regressions use inverse target-exposure variance weights.
  fit2_weights <- 1 / data_mat$se_gamma_ot^2
  fit2 <- summary(stats::lm(data_mat$Gamma_ot ~ data_mat$gamma_tr + 0,
                     weights = fit2_weights))

  ## Fit 1: weighted regression for (other) exposure association on instrument strength
  ## weights = 1 / Var(gamma_ot)
  fit1_weights <- 1 / data_mat$se_gamma_ot^2
  fit1 <- summary(stats::lm(data_mat$gamma_ot ~ data_mat$gamma_tr + 0,
                     weights = fit1_weights))

  ## Ratio-of-slopes Wald estimate
  pe <- fit2$coefficients[1] / fit1$coefficients[1]

  ## Source slope-based SE: outcome slope SE divided by exposure slope
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
#' @examples
#' x <- seq(0.08, 0.24, length.out = 24)
#' dat <- data.frame(gamma_tr = x, gamma_ot = 1.3 * x,
#'   Gamma_ot = 0.5 * 1.3 * x + 0.008 * sin(seq_along(x)),
#'   se_gamma_ot = 0.03, se_gamma_tr = 0.02, se_Gamma_ot = 0.04)
#' set.seed(2026)
#' mr_wald_bs(dat, repit = 30)
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
  sdCI   = stats::sd(result)
  return(list(pe= pe,lb =pe -sdCI*1.96 ,ub = pe + sdCI*1.96 ))

}


#' @name mr_wald_R
#' @title Median-score MR Wald estimator on a fixed grid
#' @description
#' Minimizes the squared median score over a grid with spacing 0.001 and
#' returns the source implementation's threshold-crossing interval endpoints.
#'
#' @param data_mat A data.frame containing \code{Gamma_ot}, \code{gamma_ot},
#'   \code{gamma_tr}, and \code{se_gamma_tr}. Rows correspond to harmonized SNPs;
#'   \code{gamma_tr} is the external exposure association, \code{gamma_ot} is
#'   the target exposure association, and \code{Gamma_ot} is the target outcome
#'   association.
#' @param min_num Numeric. Lower endpoint of the search grid. Default is -5.
#' @param max_num Numeric. Upper endpoint of the search grid. Default is 5.
#'
#' @details The numerical algorithm is retained from the supplied source.
#' The returned bounds are not a general implementation of the complete
#' score-inversion confidence set: disconnected accepted sets cannot be
#' represented, and an endpoint may be \code{NA} if its threshold is not
#' crossed within the grid. There is no dedicated fallback for zero score
#' variance. These limitations should be resolved before using this function
#' as a general implementation of the paper's confidence-set procedure.
#'
#' @return A named numeric vector with elements \code{pe} (a grid minimizer of
#'   the squared score), \code{lb}, and \code{ub}.
#'
#' @examples
#' x <- seq(0.08, 0.24, length.out = 24)
#' dat <- data.frame(gamma_tr = x, gamma_ot = 1.3 * x,
#'   Gamma_ot = 0.5 * 1.3 * x + 0.008 * sin(seq_along(x)),
#'   se_gamma_ot = 0.03, se_gamma_tr = 0.02, se_Gamma_ot = 0.04)
#' mr_wald_R(dat, min_num = 0.3, max_num = 0.7)
#'
#' @export
mr_wald_R <- function(data_mat,min_num = -5,max_num = 5){

  g_beta <- function(beta){
    (sum(1/data_mat$se_gamma_tr*data_mat$gamma_tr*(0.5-(data_mat$Gamma_ot > beta*data_mat$gamma_ot))))
  }
  
  beta = seq(min_num,max_num,0.001)
  M_square = unlist(lapply(beta,g_beta))^2
  index = which.min(M_square)
  pe = beta[index]
  
  CI_can  = unlist(lapply(beta,g_beta))/sqrt(0.25*(sum((data_mat$gamma_tr/data_mat$se_gamma_tr)^2)))
  CI_up   = beta[min(which(CI_can>1.96))]
  CI_low  = beta[max(which(CI_can< -1.96))]
  return(c(pe = pe,lb = CI_low,ub = CI_up))
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
#' @param alpha_star Numeric scalar or vector of length \code{p}. Added to the
#'   randomly generated direct effects. Default is 0.
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
#' @examples
#' sim <- data_gen(seed = 1, n = 200, p = 4,
#'   gamma = rep(0.1, 4), gamma_fun = function(x) 1.2 * x, beta_0 = 0.5)
#' head(sim$mat_all)
#'
#' @export
data_gen <- function(seed       = NULL,
                     n          = 10000,
                     p          = 200,
                     mu         = 0,
                     alpha_star = 0,
                     tau0 = 0,
                     gamma,
                     gamma_fun,
                     MAF        = 0.3,
                     beta_0){

  ## A supplied seed index is multiplied by 2025, as in the source.
  if (is.null(seed))
    seed <- floor(stats::runif(1) * 10000)
  set.seed(seed * 2025)

  ## Outcome sample genotype matrix: n x p, genotype in {0,1,2}
  Z <- matrix(stats::rbinom(n * p, 2, MAF), ncol = p)

  ## Unobserved confounder
  U <- stats::rnorm(n)

  ## Exposure in outcome sample
  ## NOTE: gamma_fun(gamma) must exist; otherwise replace with gamma
  D <- Z %*% gamma_fun(gamma) + U + stats::rnorm(n)

  ## Outcome in outcome sample
  ##
  alpha = stats::rnorm(p,mean  = mu , sd = tau0)+alpha_star

  Y <- beta_0 * D + U + stats::rnorm(n) + Z %*% alpha

  ## Treatment sample genotypes and exposure
  Z_new <- matrix(stats::rbinom(n * p, 2, MAF), ncol = p)
  U <- stats::rnorm(n)
  D_new <- Z_new %*% gamma + U + stats::rnorm(n)

  ## Helper: SNP -> (estimate, SE) for association with D (outcome sample)
  get_gamma_ot <- function(V){
    summary(stats::lm(D ~ V))$coefficients[2, 1:2]
  }
  gamma_ot <- apply(Z, 2, get_gamma_ot)
  se_gamma_ot <- gamma_ot[2, ]
  gamma_ot <- gamma_ot[1, ]

  ## Helper: SNP -> (estimate, SE) for association with Y (outcome sample)
  get_Gamma_ot <- function(V){
    summary(stats::lm(Y ~ V))$coefficients[2, 1:2]
  }
  Gamma_ot <- apply(Z, 2, get_Gamma_ot)
  se_Gamma_ot <- Gamma_ot[2, ]
  Gamma_ot <- Gamma_ot[1, ]

  ## Helper: SNP -> (estimate, SE) for association with D_new (treatment sample)
  get_gamma_tr <- function(V){
    summary(stats::lm(D_new ~ V))$coefficients[2, 1:2]
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


