rm(list= ls())
library(TwoSampleMR)
library(MendelianRandomization)

Sys.setenv(OPENGWAS_JWT="Use your own token after configuration")

iv_screening = extract_instruments("ieu-a-835")
# use giant dataset as screening data set that select valid IV 

IV_BMI_Southasian              = extract_outcome_data(snps = iv_screening$SNP,"ukb-e-23104_CSA",proxies = FALSE)
IV_BMI_Southasian$beta.outcome = IV_BMI_Southasian$beta.outcome*(-1)
# South Asian, the data was preprocessed using -1 transformation
# the (-1) multiplication is used so that it is comparable with other data set now
IV_BMI_Japan         = extract_outcome_data(snps = iv_screening$SNP,"bbj-a-1",proxies = FALSE)
#Japanese
IV_BMI_LatinAmerican = extract_outcome_data(snps = iv_screening$SNP,"ebi-a-GCST90095034",proxies = FALSE)
#Latin American

# can use the following code for quick check 
# plot_data = harmonise_data(convert_outcome_to_exposure(IV_BMI_Southasian),IV_BMI_Japan)
# plot(plot_data$beta.outcome~plot_data$beta.exposure)
# abline(0, 1)
# now you see all dots are clustered close to y=x
# without the -1 transformation, it will be surrounded around y=-x, which will distort some of the results


# Summary statistics for African population
# Need to download the file first !!!!! 
# <<<GCST90475155.tsv>>>>
# The file is available at https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475155/
# The summary statistics information is available at GWAS catalog: https://www.ebi.ac.uk/gwas/studies/GCST90475155

African             <- data.table::fread("*** fill in the root file here***")
afc_data = African[na.omit( match(iv_screening$SNP, (African$rsid))),]
afc_data = afc_data[,c(1:9,11)]
names(afc_data)     = c("chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","SNP","samplesize.outcome")
afc_data$outcome    = "BMI｜African"
afc_data$id.outcome = "Study: GCST90475155"
afc_data = data.frame(afc_data)
IV_BMI_African = afc_data

IV_BMI_UKB = extract_outcome_data(snps = iv_screening$SNP, "ebi-a-GCST90013974",proxies = FALSE)
IV_HDL_UKB = extract_outcome_data(snps =  iv_screening$SNP,"ebi-a-GCST90014007",proxies = FALSE)

Data1 <- list(IV_BMI_Southasian    = convert_outcome_to_exposure(IV_BMI_Southasian),
              IV_BMI_Japan         = convert_outcome_to_exposure(IV_BMI_Japan),
              IV_BMI_LatinAmerican = convert_outcome_to_exposure(IV_BMI_LatinAmerican),
              IV_BMI_African       = convert_outcome_to_exposure(IV_BMI_African),
              
              IV_BMI_UKB           = IV_BMI_UKB,
              IV_HDL_UKB           = IV_HDL_UKB) 

saveRDS(Data1,file = "Data_threshold1.rds")



iv_screening = extract_instruments("ieu-a-835",p1 = 1e-4)
# use giant dataset as screening dataset to select valid IV 

IV_BMI_Southasian    = extract_outcome_data(snps = iv_screening$SNP,"ukb-e-23104_CSA",proxies = FALSE)
#South Asian
IV_BMI_Southasian$beta.outcome = IV_BMI_Southasian$beta.outcome*(-1)


IV_BMI_Japan         = extract_outcome_data(snps = iv_screening$SNP,"bbj-a-1",proxies = FALSE)
#Japanese
IV_BMI_LatinAmerican = extract_outcome_data(snps = iv_screening$SNP,"ebi-a-GCST90095034",proxies = FALSE)
#Latin American


# Summary statistics for African population
# Need to download the file first !!!!! 
# <<<GCST90475155.tsv>>>>
# the file is available at https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475155/
# The summary statistics information is available at GWAS catalog: https://www.ebi.ac.uk/gwas/studies/GCST90475155

African             <- data.table::fread("*** fill in the root file here***")
afc_data = African[na.omit( match(iv_screening$SNP, (African$rsid))),]
afc_data = afc_data[,c(1:9,11)]
names(afc_data)     = c("chr","pos","effect_allele.outcome","other_allele.outcome","beta.outcome","se.outcome","eaf.outcome","pval.outcome","SNP","samplesize.outcome")
afc_data$outcome    = "BMI｜African"
afc_data$id.outcome = "Study: GCST90475155"
afc_data = data.frame(afc_data)
IV_BMI_African = afc_data

IV_BMI_UKB = extract_outcome_data(snps = iv_screening$SNP, "ebi-a-GCST90013974")
IV_HDL_UKB = extract_outcome_data(snps = iv_screening$SNP, "ebi-a-GCST90014007")

Data2 <- list(IV_BMI_Southasian    = convert_outcome_to_exposure(IV_BMI_Southasian),
              IV_BMI_Japan         = convert_outcome_to_exposure(IV_BMI_Japan),
              IV_BMI_LatinAmerican = convert_outcome_to_exposure(IV_BMI_LatinAmerican),
              IV_BMI_African       = convert_outcome_to_exposure(IV_BMI_African),
              IV_BMI_UKB           = IV_BMI_UKB,
              IV_HDL_UKB           = IV_HDL_UKB) 

saveRDS(Data2,file = "Data_threshold2.rds")
