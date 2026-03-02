#' @import quantreg
#' @import MASS
#' @import TwoSampleMR
#' @import mr.raps

## ----------------------------------------------------------------------------
## MR Wald-type estimator based on two weighted regressions
##   - Regress outcome association (Gamma_ot) on gamma_tr (no intercept)
##   - Regress gamma_ot on gamma_tr (no intercept)
##   - Point estimate is ratio of slopes
##   - SE uses (approx) delta-style: se(b_out) / b_exposure (your current form)
## Input:
##   data_mat must contain: Gamma_ot, gamma_ot, gamma_tr, se_Gamma_ot, se_gamma_ot
## Output:
##   pe: causal estimate; lb/ub: normal CI; bhat: fitted slope for Gamma regression

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


## ----------------------------------------------------------------------------
## Bootstrap CI for mr_wald (pairs bootstrap over SNPs)
## Input:
##   repit: number of bootstrap replications
## Output:
##   pe: original estimate; lb/ub: percentile bootstrap CI

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


## ----------------------------------------------------------------------------
## Quantile-regression analogue (median regression) for Wald-type ratio
## Note: rq() weights are treated as case weights; you used 1/se (not 1/se^2).
## Input:
##   data_mat must contain: Gamma_ot, gamma_ot, gamma_tr, se_Gamma_ot, se_gamma_ot
## Output:
##   pe: ratio of median-regression slopes

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


## ----------------------------------------------------------------------------
## Bootstrap CI for mr_wald_qr

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


## ----------------------------------------------------------------------------
## Data generator (two-sample style): generate Z,U,D,Y in an "outcome sample"
## and independently generate Z_new,U,D_new in a "treatment sample", then compute
## SNP-wise marginal associations and SEs via simple linear regressions.
##
## IMPORTANT potential issues in your current code (not changing, just flagging):
##   1) D = Z %*% g_gamma(gamma) uses g_gamma(), but you pass gamma/gamma_fun.
##      If g_gamma() is not defined elsewhere, this will error.
##   2) Y includes Z %*% alpha, but alpha is not defined in the function inputs.
##      You have alpha_star input but never used; likely you meant alpha_star.
##   3) D_new = Z_new %*% gamma uses gamma as a vector (consistent), but D uses g_gamma().
##   4) You create mat_h and mat_all with names beta.outcome/beta.exposure etc,
##      but your mr_wald() expects Gamma_ot/gamma_ot/gamma_tr column names.
##
## Output:
##   mat_h: (Gamma_outcome, gamma_tr, se_Gamma_outcome, se_gamma_tr)
##   mat_all: (gamma_ot, gamma_tr, se_gamma_ot, se_gamma_tr)

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


