rm(list= ls())
library(TwoSampleMR)
library(ieugwasr)
library(mr.raps)
library(latex2exp)
library(ggplot2)
library(dplyr)
library(quantreg)
library(MendelianRandomization)

get_number_from_string<- function(str){
  as.numeric(regmatches(str , gregexpr("-?\\d+\\.\\d+",str))[[1]])
}


getplot     <- function(str_list , estimator, population_list,clip = NULL){
  data_table = do.call(rbind,lapply(str_list,get_number_from_string))
  if(is.null(clip)){
    clip = range(data_table)
  }
  n_population = length(population_list)/2
  base_data =
    tibble::tibble(mean    = data_table[,1],
                   lower   = data_table[,2],
                   upper   = data_table[,3],
                   name    = population_list,
                   threshold = c(rep("5e-8",n_population), rep(1e-4,n_population)))
  base_data <- base_data %>%
    mutate(
      se = (upper - lower) / (1.96 * 2)
    )
  
  fig <- ggforestplot::forestplot(
    df = base_data,
    estimate = mean,
    se = se,
    colour = threshold,
    title = estimator,
    xlab = expression(hat(beta)[0]),
    xlim = range(data_table)#xlim_set
  ) +   
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  fig
}

getplotlist       <-   function(str_list_mat , clip = NULL ){
  result_list = list()
  
  xlim_set = t(sapply(str_list_mat[,1], function(x) {
    as.numeric(unlist(regmatches(x, regexec("([0-9.]+)\\(([0-9.]+),([0-9.]+)\\)", x)))[2:4])
  }))
  
  xlim_set = max(abs(xlim_set)[,2:3])*1.2
  xlim_set = c(-xlim_set,xlim_set)
  
  nfigure = ncol(str_list_mat)-1
  
  population_list = str_list_mat[,1]
  str_list_mat = str_list_mat[,-1]
  for (i in 1:nfigure) {
    result_list[[i]] <- getplot(str_list_mat[,i],names(str_list_mat)[i],population_list,clip = xlim_set)
  }
  
  result_list
}

getplot_combined  <- function(str_list_mat,Title,x_lim = NULL){
  tr_source   = str_list_mat[,1]
  n_estimator = ncol(str_list_mat)-1
  base_data = tibble::tibble()
  for (i in 1:n_estimator) {
    data_table = do.call(rbind,lapply(str_list_mat[,i+1],get_number_from_string))
    
    base_data = rbind(base_data,
                      tibble::tibble(mean      = data_table[,1],
                                     lower     = data_table[,2],
                                     upper     = data_table[,3],
                                     name      = colnames(str_list_mat)[i+1],
                                     source    = tr_source)
    )
  }
  base_data <- base_data %>%
    mutate(
      se = (upper - lower) / (1.96 * 2)
    )
  
  if(is.null(x_lim))
    x_lim = range(range(base_data[,1:3]))
  
  fig <- ggforestplot::forestplot(
    df = base_data,
    estimate = mean,
    se = se,
    colour = source,
    title = Title,
    xlab = expression(hat(beta)[0]),
    xlim = x_lim#xlim_set
  ) +   
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    labs(color = "Treatment Data Source")
  
  
  fig
  
  
}

getresultandci_lb_ub <- function(pe,lb,ub){
  vec  = c(pe,lb,ub)
  vec = round(vec,digit = 3)
  paste0(vec[1],"(",vec[2],",",vec[3],")")
}

mr_ivqr        <- function(data_mat){
  g_beta <- function(beta){
    (sum(1/data_mat$se_gamma_tr*data_mat$gamma_tr*(0.5-(data_mat$Gamma_ot > beta*data_mat$gamma_ot))))
  }
  
  min_num = range(data_mat$Gamma_ot/data_mat$gamma_ot)[1]
  max_num = range(data_mat$Gamma_ot/data_mat$gamma_ot)[2]
  
  beta = seq(min_num,max_num,0.001)
  M_square = unlist(lapply(beta,g_beta))^2
  index = which.min(M_square)
  pe = beta[index]
  
  CI_can  = unlist(lapply(beta,g_beta))/sqrt(0.25*( sum(1+(data_mat$gamma_tr/data_mat$se_gamma_tr)^2)))
  CI_up   = beta[min(which(CI_can>1.96))]
  CI_low  = beta[max(which(CI_can< -1.96))]
  return(c(pe = pe,lb = CI_low,ub = CI_up))
}


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
  wald_R_fit       <- RobustMR::mr_wald_R(data_full)
  
  
  final_result = data.frame(MR_Egger = c(mr.fit$b[1],mr.fit$b[1] - 1.96*mr.fit$se[1],mr.fit$b[1] + 1.96*mr.fit$se[1]))
  
  final_result$W_Median =  c(mr.fit$b[2],mr.fit$b[2] - 1.96*mr.fit$se[2],mr.fit$b[2] + 1.96*mr.fit$se[2])
  final_result$IVW      =  c(mr.fit$b[3],mr.fit$b[3] - 1.96*mr.fit$se[3],mr.fit$b[3] + 1.96*mr.fit$se[3])
  
  final_result$DIVW     = c(divwfit@Estimate, divwfit@Estimate - 1.96*divwfit@StdError,  divwfit@Estimate + 1.96*divwfit@StdError)
  final_result$RAPS     = c(mr.raps.fit$beta.hat, mr.raps.fit$beta.hat - 1.96*mr.raps.fit$beta.se, mr.raps.fit$beta.hat + 1.96*mr.raps.fit$beta.se)
  
  final_result$MR_Wald   = c(wald_fit$pe,wald_fit$lb,wald_fit$ub)
  final_result$MR_Wald_R = wald_R_fit
  final_result
  
  
}

setwd("***use the actual address after download the data from the package*****")
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


result_mat1 = matrix(NA,ncol = 7,nrow =4 )
for (i in 1:4) {
  for (j in 1:length(result1[[1]])) {
    result_mat1[i,j] = getresultandci_lb_ub(result1[[i]][[j]][1],result1[[i]][[j]][2],result1[[i]][[j]][3])
  }
}

colnames(result_mat1) =  c("MR-Egger",  "W.Median",  "IVW" ,      "DIVW" ,     "RAPS" ,     "MR-Wald"  ,
                           "MR-Wald-R")
rownames(result_mat1) = c("South Asian","Japanese","Latin American","African")
fit1_data = cbind(data.frame(Tr_source = c("South Asian","Japanese","Latin American","African")),result_mat1)
fig1      = getplot_combined(fit1_data,"Causal effect of BMI on HDL",x_lim = c(-2.2,0.3))



result_mat2 = matrix(NA,ncol = 7,nrow =4 )
colnames(result_mat2) =  names(result2[[1]])
for (i in 1:4) {
  for (j in 1:length(result2[[1]])) {
    result_mat2[i,j] = getresultandci_lb_ub(result2[[i]][[j]][1],result2[[i]][[j]][2],result2[[i]][[j]][3])
  }
}

rownames(result_mat2) = c("South Asian","Japanese","Latin American","African")
colnames(result_mat2) =  c("MR-Egger",  "W.Median",  "IVW" ,      "DIVW" ,     "RAPS" ,     "MR-Wald"  ,
                           "MR-Wald-R")

fit2_data = cbind(data.frame(Tr_source = c("South Asian","Japanese","Latin American","African")),result_mat2)
fig2 = getplot_combined(fit2_data,"Causal effect of BMI on HDL",x_lim = c(-2.2,0.3))


ggsave("bygroup2.png",fig2+ theme(legend.position = "none"),height = 4,width = 4,units = "in")
ggsave("bygroup1.png",fig1+ theme(legend.position = "none"),height = 4,width = 4,units = "in")
