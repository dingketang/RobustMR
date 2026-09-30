small_summary_data <- function() {
  j <- seq_len(40L)
  gamma_ot <- 0.10 + 0.005 * j
  data.frame(
    Gamma_ot = 0.5 * gamma_ot + 0.005 * cos(j),
    gamma_ot = gamma_ot,
    gamma_tr = 0.08 + 0.004 * j,
    se_Gamma_ot = rep(0.03, length(j)),
    se_gamma_tr = 0.01 + 0.0001 * j,
    se_gamma_ot = 0.02 + 0.0002 * j
  )
}
