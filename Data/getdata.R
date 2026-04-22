rm(list= ls())
library(TwoSampleMR)
library(MendelianRandomization)

Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJkdGFuZ0B1b3R0YXdhLmNhIiwiaWF0IjoxNzc2ODgwNTM3LCJleHAiOjE3NzgwOTAxMzd9.UtTkOyjpCmQSO75pZkpO2CO4-eucUO_8p_yfF2M2EV5g58zTIyYU08PFq-sxf6zjZAqni5IyTuN5h4XZrEmqFo8uUVKXmAI4sSvwB0RGB8GOWpCj78jybTFHXNOLPNOmNHtyQWU4xhNLRVxx3qiL6Pju9lUhK9a_O7iexEzy2PLCRhaQUVxypvJYswVD7TDNhs95G4QejSEPDIZIU93IT1qqvqAT2HNEI9Ap-0vZexDJkdK0lXDT8PT_WTx7wemtCZ81O5cIZO_kEfCl4nqET3XyR6RGK7KEnalQhRUA0VR5pyo80D5fK3Qz99bDRcXsLkQaKKUlCkPgLnkBsSo2Gg")

iv_screening = extract_instruments("ieu-a-835")
# use giant dataset as screening data set that select valid IV 

IV_BMI_Southasian              = extract_outcome_data(snps = iv_screening$SNP,"ukb-e-23104_CSA",proxies = FALSE)
IV_BMI_Southasian$beta.outcome = IV_BMI_Southasian$beta.outcome*(-1)
# South Asian, the data was prepossessed using -1 transformation
# the (-1) multiplication is used so that it is comparable with other data set now
IV_BMI_Japan         = extract_outcome_data(snps = iv_screening$SNP,"bbj-a-1",proxies = FALSE)
#Japanese
IV_BMI_LatinAmerican = extract_outcome_data(snps = iv_screening$SNP,"ebi-a-GCST90095034",proxies = FALSE)
#Latin American

# can use the following code for quick check 
# plot_data = harmonise_data(convert_outcome_to_exposure(IV_BMI_Southasian),IV_BMI_Japan)
# plot(plot_data$beta.outcome~plot_data$beta.exposure)
# abline(0, 1)
# now you see all dots are surrounded close to y=x
# without the -1 transfermation, it will be surranded arround y=-x, which will distort some of the result


# Summary statistics for African population
# Need to download the file first !!!!! 
# <<<GCST90475155.tsv>>>>
# the file is available at https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475155/
# the summary statistics information are available at GWAS catalog: https://www.ebi.ac.uk/gwas/studies/GCST90475155

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
# use giant dataset as screening data set that select valid IV 

IV_BMI_Southasian    = extract_outcome_data(snps = iv_screening$SNP,"ukb-e-23104_CSA",proxies = FALSE)
#South Asian
IV_BMI_Southasian$beta.outcome = IV_BMI_Southasian$beta.outcome*(-1)


IV_BMI_Japan         = extract_outcome_data(snps = iv_screening$SNP,"bbj-a-1",proxies = FALSE)
#Japanese
IV_BMI_LatinAmerican = extract_outcome_data(snps = iv_screening$SNP,"ebi-a-GCST90095034",proxies = FALSE)
#Latin American


# Summary statistics for african population
# Need to download the file first !!!!! 
# <<<GCST90475155.tsv>>>>
# the file is available at https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475155/
# the summary statistics information are available at GWAS catalog: https://www.ebi.ac.uk/gwas/studies/GCST90475155

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
