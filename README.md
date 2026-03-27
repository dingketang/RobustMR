# RobustMR

This repository provides the R code and package accompanying the paper  
**“A Robust Framework for Two-Sample Mendelian Randomization under Population Heterogeneity.”**

---

## Installation

You can install the **RobustMR** package directly from GitHub:

```r
devtools::install_github("dingketang/RobustMR")
```

---

## Dependencies

The **RobustMR** package depends on the following R packages. The versions listed below correspond to those used when building and testing the package.

- quantreg (6.1)
- TwoSampleMR (0.6.11)
- ieugwasr (1.1.0)
- mr.raps (0.4.2)
- parallel (4.4.1)
- MASS (7.3-65)

---

## Example Usage

The following example illustrates how to generate simulated data and apply the proposed MR estimators.

```r
n = 10000  # sample size
p = 200    # number of SNPs

# IV–treatment effects in the treatment GWAS
set.seed(1234)
gamma = runif(p, min = 0.05, max = 0.1)

gamma_fun = function(gamma) {
  (gamma + 0.1) / 2
}

# Heterogeneous setting with a linear transformation
tau0 = 0
mu   = 0
alpha_star = rep(0, p)  # no pleiotropy

data = RobustMR::data_gen(
  seed = 1234,
  n = n,
  p = p,
  mu = mu,
  alpha_star = alpha_star,
  tau0 = tau0,
  gamma = gamma,
  gamma_fun = gamma_fun,
  MAF = 0.3,
  beta_0 = 0.5
)

set.seed(1234)
RobustMR::mr_wald_bs(data$mat_all)
# Point estimate ≈ 0.511
# 95% CI ≈ (0.473, 0.549)

set.seed(1234)
RobustMR::mr_wald_qr_bs(data$mat_all)
# Point estimate ≈ 0.501
# 95% CI ≈ (0.445, 0.557)
```

---

## Reproducing Results from the Manuscript

To facilitate reproducibility, we provide the function `senario1case3`, which replicates part of the simulation results corresponding to:

- Scenario (i): heterogeneous setting  
- Case (3): sine transformation  

By modifying this function, the other scenarios reported in the manuscript can be reproduced.

```r
## ----------------------------------------------------------------------------
## Run Monte Carlo experiment (1000 replications, 5 cores)

result <- mclapply(1:1000, senario1case3, mc.cores = 5)

# Runtime is approximately 10 minutes on a MacBook Pro (2021)

## Summarize performance
process_fit_result(fitlist = result, beta_0 = 0.5)

## ----------------------------------------------------------------------------
#          mr_wald mr_wald_qr mr_w_median mr_egger mr_rap mr_divw
# Bias          0.0        0.0        75.9    -59.2  136.7   137.3
# RMSE          1.8        3.3        76.4     62.4  137.0   137.6
# CI length     7.0       13.8        41.1     67.1   35.0    33.7
# CI           94.4       95.2         0.0     11.2    0.0     0.0
```

---

## Notes

- Parallel computing is supported via the `parallel` package.
- Simulation results may vary slightly due to randomness and system differences.
