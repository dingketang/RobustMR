rm(list= ls())
library(TwoSampleMR)
library(ieugwasr)
library(mr.raps)
library(latex2exp)
library(ggplot2)
library(dplyr)
library(quantreg)
library(MendelianRandomization)

MR_all<- function(Gamma_ou,gamma_ou, gamma_tr){
  dat1 = harmonise_data(gamma_tr,Gamma_ou)
  dat2 = harmonise_data(gamma_tr,gamma_ou)
  
  mr.fit  = mr(dat1,method_list=c("mr_egger_regression","mr_weighted_median", "mr_ivw"))
  divwfit = MendelianRandomization::mr_divw(mr_input(bx   = dat1$beta.exposure,
                                                     by   = dat1$beta.outcome,
                                                     bxse = dat1$se.exposure,
                                                     byse = dat1$se.outcome))
  
  mr.raps.fit <- mr.raps(dat1,diagnostics = FALSE)
  
  intersect_vals <- intersect(dat1$SNP, dat2$SNP)
  
  indices_vec1 <- which(dat1$SNP %in% intersect_vals)
  indices_vec2 <- which(dat2$SNP %in% intersect_vals)
  
  dat1_new <- dat1[indices_vec1,]
  dat2_new <- dat2[indices_vec2,]
  
  data_full = data.frame(gamma_tr    = dat1_new$beta.exposure,
                         se_gamma_tr = dat1_new$se.exposure,
                         Gamma_ot    = dat1_new$beta.outcome,
                         se_Gamma_ot = dat1_new$se.outcome,
                         gamma_ot    = dat2_new$beta.outcome,
                         se_gamma_ot = dat2_new$se.outcome)
  
  wald_fit         <- RobustMR::mr_wald_bs(data_full)
  wald_qr_fit      <- RobustMR::mr_wald_qr_bs(data_full)
  
  final_result = data.frame(MR_Egger = c(mr.fit$b[1],mr.fit$b[1] - 1.96*mr.fit$se[1],mr.fit$b[1] + 1.96*mr.fit$se[1]))
  
  final_result$W_Median =  c(mr.fit$b[2],mr.fit$b[2] - 1.96*mr.fit$se[2],mr.fit$b[2] + 1.96*mr.fit$se[2])
  final_result$IVW      =  c(mr.fit$b[3],mr.fit$b[3] - 1.96*mr.fit$se[3],mr.fit$b[3] + 1.96*mr.fit$se[3])
  
  final_result$DIVW     = c(divwfit@Estimate, divwfit@Estimate - 1.96*divwfit@StdError,  divwfit@Estimate + 1.96*divwfit@StdError)
  final_result$RAPS     = c(mr.raps.fit$beta.hat, mr.raps.fit$beta.hat - 1.96*mr.raps.fit$beta.se, mr.raps.fit$beta.hat + 1.96*mr.raps.fit$beta.se)
  
  final_result$MR_Wald   = c(wald_fit$pe,wald_fit$lb,wald_fit$ub)
  final_result$MR_Wald_R = c(wald_qr_fit$pe,wald_qr_fit$lb,wald_qr_fit$ub)
  
  final_result
  
  
}


data1 = readRDS("Data_threshold1.rds")
data2 = readRDS("Data_threshold2.rds")  

Tr_data_names <- names(data1)[1:4]

result1 = list()
result2 = list()

for (a in Tr_data_names){
  Gamma_ou = data1$IV_HDL_UKB
  gamma_ou = data1$IV_BMI_UKB
  gamma_tr = data1[[a]]
  fitresult <- MR_all(Gamma_ou,gamma_ou,gamma_tr)
  rownames(fitresult) = c("pe","lb","up")
  result1[[a]] = fitresult
  
  
  Gamma_ou = data2$IV_HDL_UKB
  gamma_ou = data2$IV_BMI_UKB
  gamma_tr = data2[[a]]
  fitresult <- MR_all(Gamma_ou,gamma_ou,gamma_tr)
  rownames(fitresult) = c("pe","lb","up")
  result2[[a]] = fitresult
}