# Run after installing the package.
library(RobustMR)

x <- seq(0.08, 0.24, length.out = 24)
dat <- data.frame(
  gamma_tr = x,
  gamma_ot = 1.3 * x,
  Gamma_ot = 0.5 * 1.3 * x + 0.008 * sin(seq_along(x)),
  se_gamma_tr = 0.02,
  se_gamma_ot = 0.03,
  se_Gamma_ot = 0.04
)

print(mr_wald(dat))

set.seed(2026)
print(mr_wald_bs(dat, repit = 100))

print(mr_wald_R(dat, min_num = 0.3, max_num = 0.7))
