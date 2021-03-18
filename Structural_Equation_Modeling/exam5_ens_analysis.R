##############################################################################
## Title: IL1RN Models with normalized gene expression variables
## Version: 1
## Author: Regina Manansala
## Date Created: 17-June-2020
## Date Modified: 04-January-2021
## Results File: exam1_ens analysis.docx; exam5_ens_analysis.docx
## Output Data: exam1_MESA_norm.docx; exam5_MESA_norm.docx; 
##				exam1_norm_wENS.txt; exam5_norm_wENS.txt
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

# vars <- c("IL.1.R.AcP", "IL.6", "IL.1.sRI", "IL.1b", "IL.1Ra", "IL.1a")
# summary(exam1_MESA[, vars])
# vars_norm <- c("IL.1.R.AcP_norm", "IL.6_norm", "IL.1.sRI_norm", "IL.1b_norm", "IL.1Ra_norm", "IL.1a_norm")
# summary(exam1_MESA[, vars_norm])
# for(i in vars_norm){
#   hist(exam1_MESA[[i]], main = paste0("Histogram of ", i), xlab = paste0("exam1_MESA$", i))
#   boxplot(exam1_MESA[[i]], main = paste0("Boxplot of ", i), xlab = paste0("exam1_MESA$", i))
# }
# 
# summary(lm(IL.6_norm ~ chr2.113083453_G, data=exam1_MESA))
# summary(lm(IL.1Ra_norm ~ chr2.113083453_G, data=exam1_MESA))
# summary(lm(IL.1b_norm ~ chr2.113083453_G, data=exam1_MESA))
# summary(lm(IL.1.sRI_norm ~ chr2.113083453_G, data=exam1_MESA))
# summary(lm(IL.1.R.AcP_norm ~ chr2.113083453_G, data=exam1_MESA))
# summary(lm(IL.1Ra_norm ~ IL.1b_norm, data=exam1_MESA))
# summary(lm(IL.1Ra_norm ~ IL.1.sRI_norm, data=exam1_MESA))
# summary(lm(IL.1Ra_norm ~ IL.1.R.AcP_norm, data=exam1_MESA))
# summary(lm(IL.6_norm ~ IL.1b_norm, data=exam1_MESA))
# summary(lm(IL.6_norm ~ IL.1.sRI_norm, data=exam1_MESA))
# summary(lm(IL.6_norm ~ IL.1.R.AcP_norm, data=exam1_MESA))
# summary(lm(IL.6_norm ~ IL.1Ra_norm, data=exam1_MESA))

summary(lm(IL.1a_norm ~ chr2.113083453_G, data=exam1_MESA))
summary(lm(IL.6_norm ~ IL.1a_norm, data=exam1_MESA))

#SNP vs ENS variables in Exam 1
summary(lm(ENS_PBMC_norm ~ chr2.113083453_G, data=exam1_MESA))

#ENS variables vs IL1Ra_norm in Exam 1K
summary(lm(IL.1a_norm ~ ENS_PBMC_norm, data=exam1_MESA))

#ENS variables vs IL6_norm in Exam 1
summary(lm(IL.6_norm ~ ENS_PBMC_norm, data=exam1_MESA))


# summary(lm(chr2.113083453_G ~ IL.1.R.AcP, data=exam1_MESA))
# summary(lm(chr2.113083453_G ~ IL.6, data=exam1_MESA))
# summary(lm(chr2.113083453_G ~ IL.1.sRI, data=exam1_MESA))
# summary(lm(chr2.113083453_G ~ IL.1b, data=exam1_MESA))
# summary(lm(chr2.113083453_G ~ IL.1Ra, data=exam1_MESA))

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

#SNP vs ENS variables in Exam 5
summary(lm(ENS_PBMC_norm ~ chr2.113083453_G, data=exam5_MESA))
summary(lm(ENS_MONO_norm ~ chr2.113083453_G, data=exam5_MESA))
summary(lm(ENS_TCELL_norm ~ chr2.113083453_G, data=exam5_MESA))

#ENS variables vs IL1Ra_norm in Exam 5K
summary(lm(IL.1a_norm ~ ENS_PBMC_norm, data=exam5_MESA))
summary(lm(IL.1a_norm ~ ENS_MONO_norm, data=exam5_MESA))
summary(lm(IL.1a_norm ~ ENS_TCELL_norm, data=exam5_MESA))

#ENS variables vs IL6_norm in Exam 5
summary(lm(IL.6_norm ~ ENS_PBMC_norm, data=exam5_MESA))
summary(lm(IL.6_norm ~ ENS_MONO_norm, data=exam5_MESA))
summary(lm(IL.6_norm ~ ENS_TCELL_norm, data=exam5_MESA))

# summary(exam5_MESA[, vars])
# summary(exam5_MESA[, vars_norm])
# for(i in vars_norm){
#   hist(exam5_MESA[[i]], main = paste0("Histogram of ", i), xlab = paste0("exam5_MESA$", i))
#   boxplot(exam5_MESA[[i]], main = paste0("Boxplot of ", i), xlab = paste0("exam5_MESA$", i))
# }
# 
# summary(lm(IL.6_norm ~ chr2.113083453_G, data=exam5_MESA))
# summary(lm(IL.1Ra_norm ~ chr2.113083453_G, data=exam5_MESA))
# summary(lm(IL.1b_norm ~ chr2.113083453_G, data=exam5_MESA))
# summary(lm(IL.1.sRI_norm ~ chr2.113083453_G, data=exam5_MESA))
# summary(lm(IL.1.R.AcP_norm ~ chr2.113083453_G, data=exam5_MESA))
# summary(lm(IL.1Ra_norm ~ IL.1b_norm, data=exam5_MESA))
# summary(lm(IL.1Ra_norm ~ IL.1.sRI_norm, data=exam5_MESA))
# summary(lm(IL.1Ra_norm ~ IL.1.R.AcP_norm, data=exam5_MESA))
# summary(lm(IL.6_norm ~ IL.1b_norm, data=exam5_MESA))
# summary(lm(IL.6_norm ~ IL.1.sRI_norm, data=exam5_MESA))
# summary(lm(IL.6_norm ~ IL.1.R.AcP_norm, data=exam5_MESA))
# summary(lm(IL.6_norm ~ IL.1Ra_norm, data=exam5_MESA))

summary(lm(IL.1a_norm ~ chr2.113083453_G, data=exam5_MESA))
summary(lm(IL.6_norm ~ IL.1a_norm, data=exam5_MESA))
# 
# summary(lm(chr2.113083453_G ~ IL.1.R.AcP, data=exam5_MESA))
# summary(lm(chr2.113083453_G ~ IL.6, data=exam5_MESA))
# summary(lm(chr2.113083453_G ~ IL.1.sRI, data=exam5_MESA))
# summary(lm(chr2.113083453_G ~ IL.1b, data=exam5_MESA))
# summary(lm(chr2.113083453_G ~ IL.1Ra, data=exam5_MESA))
# 
# CHS_lv <- CHS_addtl %>% mutate(., unique_subject_key = paste0("CHS_", Individual_ID)) %>% 
#   inner_join(., samples[, c("unique_subject_key", "sample.id")], by = "unique_subject_key") %>%
#   inner_join(., rs6734238[, c("IID", "chr2:113083453_G")], by = c("sample.id" = "IID")) %>%
#   inner_join(., IL6_CHS_encore, by = "sample.id") %>%
#   inner_join(., crp_encore[crp_encore$study.x == "CHS", c("sample.id", "crp")], by = "sample.id") %>%
#   mutate(., il1ray2 = log(il1ray2), il1ray5 = log(il1ray5)) %>%
#   rename(., IID = sample.id, chr2.113083453_G = `chr2:113083453_G`, IL.1Ra = il1ray2)
# CHS_lv2 <- subset(CHS_lv, CHS_lv$il6harmonizedy2 != 0 & CHS_lv$il6harmonizedy2 < mean(CHS_lv$il6harmonizedy2, na.rm = TRUE) + 3*sd(CHS_lv$il6harmonizedy2, na.rm = TRUE)) %>%
#   mutate(., il6harmonizedy2 = log(il6harmonizedy2))

## Original IL6 Variable
# summary(lm(chr2.113083453_G ~ il6, data=CHS_lv))

## IL6Harmonized2 - after removing outliers and missings
# summary(lm(chr2.113083453_G ~ il6harmonizedy2, data=CHS_lv2))


# write.table(exam1_MESA, "exam1_MESA_norm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(exam5_MESA, "exam5_MESA_norm.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


