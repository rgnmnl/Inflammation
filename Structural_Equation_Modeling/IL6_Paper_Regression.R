##############################################################################
## Title: Models for IL6 Paper
## Version: 1
## Author: Regina Manansala
## Date Created: 28-December-2020
## Date Modified: 28-December-2020
## Results File: IL6_Paper_Regression.docx
##############################################################################

library(tidyr)
library(ggplot2)
library(qqman)
library(data.table)
library(dplyr)
library(lavaan)

data_dir <- "~/Documents/Inflammation_SEM/IL6_Analysis/"
setwd(data_dir)

rs6734238 <- fread("freeze8_rs6734238.raw")
exam1 <- fread("IL1RN_phenotypes_exam1.txt")
exam5 <- fread("IL1RN_phenotypes_exam5.txt")
CHS_addtl <- fread("CHS_MLB_20200121_CHIP_INFLAMMATION.csv")
samples <- fread("~/Documents/Inflammation_Data/freeze8_sample_annot_2019-03-28.txt")

## Make Mediation Analysis Data

# crp_encore <- fread("~/Documents/Inflammation_Data_Results/CRP/log_crp_encore.txt")
crp_encore <- fread("~/Documents/Inflammation_Data_Results/CRP/log_crp_encore.txt")
crp_encore$study.x <- gsub("_.*", "", crp_encore$race_study)
IL6_MESA_encore <- fread("~/Documents/Inflammation_Data_Results/IL6/il6_MESA_20190826.txt")
IL6_CHS_encore <- fread("~/Documents/Inflammation_Data_Results/IL6/il6_CHS_20190826.txt")

################ Make Exam 1 data set ################ 

exam1_MESA <- exam1 %>% inner_join(., IL6_MESA_encore, by = c("IID" = "sample.id")) %>% 
  inner_join(., crp_encore[crp_encore$study.x == "MESA", c("sample.id", "crp")], by = c("IID" = "sample.id")) %>%
  mutate(., IL.1Ra = log(IL.1Ra), IL.1.R.AcP = log(IL.1.R.AcP), IL.6 = log(IL.6), IL.1.sRI = log(IL.1.sRI), IL.1b = log(IL.1b), IL.1a = log(IL.1a))

# Delete observations with log(IL6) > 8
exam1_MESA$IL.6 <- ifelse(exam1_MESA$IL.6 > 7, NA, exam1_MESA$IL.6)
# Delete observations with log(IL1B) > 10
exam1_MESA$IL.1b <- ifelse(exam1_MESA$IL.1b > 10, NA, exam1_MESA$IL.1b)
# delete log(IL1Ra) > 10
exam1_MESA$IL.1Ra <- ifelse(exam1_MESA$IL.1Ra > 10, NA, exam1_MESA$IL.1Ra)
# Delete observations with log(IL6) > 8
exam1_MESA$IL.1a <- ifelse(exam1_MESA$IL.1a > 8, NA, exam1_MESA$IL.1a)

# normalize the values (subtract the mean and divide by the standard deviation).
exam1_MESA$IL.1.R.AcP_norm <- (exam1_MESA$IL.1.R.AcP - mean(exam1_MESA$IL.1.R.AcP, na.rm = TRUE)) / sd(exam1_MESA$IL.1.R.AcP, na.rm = TRUE)
exam1_MESA$IL.6_norm <- (exam1_MESA$IL.6 - mean(exam1_MESA$IL.6, na.rm = TRUE)) / sd(exam1_MESA$IL.6, na.rm = TRUE)
exam1_MESA$IL.1.sRI_norm <- (exam1_MESA$IL.1.sRI - mean(exam1_MESA$IL.1.sRI, na.rm = TRUE)) / sd(exam1_MESA$IL.1.sRI, na.rm = TRUE)
exam1_MESA$IL.1b_norm <- (exam1_MESA$IL.1b - mean(exam1_MESA$IL.1b, na.rm = TRUE)) / sd(exam1_MESA$IL.1b, na.rm = TRUE)
exam1_MESA$IL.1Ra_norm <- (exam1_MESA$IL.1Ra - mean(exam1_MESA$IL.1Ra, na.rm = TRUE)) / sd(exam1_MESA$IL.1Ra, na.rm = TRUE)
exam1_MESA$IL.1a_norm <- (exam1_MESA$IL.1a - mean(exam1_MESA$IL.1a, na.rm = TRUE)) / sd(exam1_MESA$IL.1a, na.rm = TRUE)

#Normalize ENS variables
exam1_MESA$ENSG00000136689.18_PBMC_exam1 <- log(exam1_MESA$ENSG00000136689.18_PBMC_exam1)
exam1_MESA$ENS_PBMC_norm <- (exam1_MESA$ENSG00000136689.18_PBMC_exam1 - mean(exam1_MESA$ENSG00000136689.18_PBMC_exam1, na.rm = TRUE)) / sd(exam1_MESA$ENSG00000136689.18_PBMC_exam1, na.rm = TRUE)

################ Make Exam 5 data set ################ 

exam5_MESA <- exam5 %>% inner_join(., IL6_MESA_encore, by = c("IID" = "sample.id")) %>% 
  inner_join(., crp_encore[crp_encore$study.x == "MESA", c("sample.id", "crp")], by = c("IID" = "sample.id")) %>%
  mutate(., IL.1Ra = log(IL.1Ra), IL.1.R.AcP = log(IL.1.R.AcP), IL.6 = log(IL.6), IL.1.sRI = log(IL.1.sRI), IL.1b = log(IL.1b), IL.1a = log(IL.1a))

# Delete observations with log(IL6) > 8
exam5_MESA$IL.6 <- ifelse(exam5_MESA$IL.6 > 7, NA, exam5_MESA$IL.6)
# Delete observations with log(IL1B) > 10
exam5_MESA$IL.1b <- ifelse(exam5_MESA$IL.1b > 10, NA, exam5_MESA$IL.1b)
# delete log(IL1Ra) > 10
exam5_MESA$IL.1Ra <- ifelse(exam5_MESA$IL.1Ra > 10, NA, exam5_MESA$IL.1Ra)
# Delete observations with log(IL6) > 8
exam5_MESA$IL.1a <- ifelse(exam5_MESA$IL.1a > 8, NA, exam5_MESA$IL.1a)

# normalize the values (subtract the mean and divide by the standard deviation).
exam5_MESA$IL.1.R.AcP_norm <- (exam5_MESA$IL.1.R.AcP - mean(exam5_MESA$IL.1.R.AcP, na.rm = TRUE)) / sd(exam5_MESA$IL.1.R.AcP, na.rm = TRUE)
exam5_MESA$IL.6_norm <- (exam5_MESA$IL.6 - mean(exam5_MESA$IL.6, na.rm = TRUE)) / sd(exam5_MESA$IL.6, na.rm = TRUE)
exam5_MESA$IL.1.sRI_norm <- (exam5_MESA$IL.1.sRI - mean(exam5_MESA$IL.1.sRI, na.rm = TRUE)) / sd(exam5_MESA$IL.1.sRI, na.rm = TRUE)
exam5_MESA$IL.1b_norm <- (exam5_MESA$IL.1b - mean(exam5_MESA$IL.1b, na.rm = TRUE)) / sd(exam5_MESA$IL.1b, na.rm = TRUE)
exam5_MESA$IL.1Ra_norm <- (exam5_MESA$IL.1Ra - mean(exam5_MESA$IL.1Ra, na.rm = TRUE)) / sd(exam5_MESA$IL.1Ra, na.rm = TRUE)
exam5_MESA$IL.1a_norm <- (exam5_MESA$IL.1a - mean(exam5_MESA$IL.1a, na.rm = TRUE)) / sd(exam5_MESA$IL.1a, na.rm = TRUE)

#Normalize ENS variables
exam5_MESA$ENSG00000136689.18_PBMC_exam5 <- log(exam5_MESA$ENSG00000136689.18_PBMC_exam5)
exam5_MESA$ENS_PBMC_norm <- (exam5_MESA$ENSG00000136689.18_PBMC_exam5 - mean(exam5_MESA$ENSG00000136689.18_PBMC_exam5, na.rm = TRUE)) / sd(exam5_MESA$ENSG00000136689.18_PBMC_exam5, na.rm = TRUE)

exam5_MESA$ENSG00000136689.18_Mono_exam5 <- log(exam5_MESA$ENSG00000136689.18_Mono_exam5)
exam5_MESA$ENS_MONO_norm <- (exam5_MESA$ENSG00000136689.18_Mono_exam5 - mean(exam5_MESA$ENSG00000136689.18_Mono_exam5, na.rm = TRUE)) / sd(exam5_MESA$ENSG00000136689.18_Mono_exam5, na.rm = TRUE)

exam5_MESA$ENSG00000136689.18_Tcell_exam5 <- ifelse(exam5_MESA$ENSG00000136689.18_Tcell_exam5 == 0, 0.0001, exam5_MESA$ENSG00000136689.18_Tcell_exam5)
exam5_MESA$ENSG00000136689.18_Tcell_exam5 <- log(exam5_MESA$ENSG00000136689.18_Tcell_exam5)
exam5_MESA$ENS_TCELL_norm <- (exam5_MESA$ENSG00000136689.18_Tcell_exam5 - mean(exam5_MESA$ENSG00000136689.18_Tcell_exam5, na.rm = TRUE)) / sd(exam5_MESA$ENSG00000136689.18_Tcell_exam5, na.rm = TRUE)

############################################################

############## EXAM 1 ##############

summary(exam1_MESA$ENS_PBMC_norm)

#SNP vs ENS variables in Exam 1
summary(lm(ENS_PBMC_norm ~ chr2.113083453_G, data=exam1_MESA))

#ENS variables vs IL1Ra_norm in Exam 1
summary(lm(IL.1Ra_norm ~ ENS_PBMC_norm, data=exam1_MESA))

#ENS variables vs IL6_norm in Exam 1
summary(lm(IL.6_norm ~ ENS_PBMC_norm, data=exam1_MESA))

#ENS vs SNP in Exam 1 - adjusted
summary(lm(ENS_PBMC_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam1_MESA))

#IL6 protein vs SNP in Exam 1 - adjusted
summary(lm(IL.6_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam1_MESA))

#IL1Ra protein vs SNP in Exam 1 - adjusted
summary(lm(IL.1Ra_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam1_MESA))

############## EXAM 5 ##############

summary(exam5_MESA$ENS_MONO_norm)
summary(exam5_MESA$ENS_PBMC_norm)
summary(exam5_MESA$ENS_TCELL_norm)

#SNP vs ENS variables in Exam 5
summary(lm(ENS_PBMC_norm ~ chr2.113083453_G, data=exam5_MESA))
summary(lm(ENS_MONO_norm ~ chr2.113083453_G, data=exam5_MESA))
summary(lm(ENS_TCELL_norm ~ chr2.113083453_G, data=exam5_MESA))

#ENS variables vs IL1Ra_norm in Exam 5K
summary(lm(IL.1Ra_norm ~ ENS_PBMC_norm, data=exam5_MESA))
summary(lm(IL.1Ra_norm ~ ENS_MONO_norm, data=exam5_MESA))
summary(lm(IL.1Ra_norm ~ ENS_TCELL_norm, data=exam5_MESA))

#ENS variables vs IL6_norm in Exam 5
summary(lm(IL.6_norm ~ ENS_PBMC_norm, data=exam5_MESA))
summary(lm(IL.6_norm ~ ENS_MONO_norm, data=exam5_MESA))
summary(lm(IL.6_norm ~ ENS_TCELL_norm, data=exam5_MESA))

#ENS variables vs SNP in Exam 5 - adjusted
summary(lm(ENS_PBMC_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam5_MESA))
summary(lm(ENS_MONO_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam5_MESA))
summary(lm(ENS_TCELL_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam5_MESA))

#IL6 protein vs SNP in Exam 5 - adjusted
summary(lm(IL.6_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam5_MESA))

#IL1Ra protein vs SNP in Exam 5 - adjusted
summary(lm(IL.1Ra_norm ~ chr2.113083453_G + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=exam5_MESA))







