##############################################################################
## Title: IL1Rn Additional Models
## Version: 1
## Author: Regina Manansala
## Date Created: 21-April-2020
## Date Modified: 22-April-2020
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
  mutate(., IL.1Ra = log(IL.1Ra), IL.1.R.AcP = log(IL.1.R.AcP), IL.6 = log(IL.6), IL.1.sRI = log(IL.1.sRI), IL.1b = log(IL.1b))

# Delete observations with log(IL6) > 8
exam1_MESA$IL.6 <- ifelse(exam1_MESA$IL.6 > 7, NA, exam1_MESA$IL.6)
# Delete observations with log(IL1B) > 10
exam1_MESA$IL.1b <- ifelse(exam1_MESA$IL.1b > 10, NA, exam1_MESA$IL.1b)
# delete log(IL1Ra) > 10
exam1_MESA$IL.1Ra <- ifelse(exam1_MESA$IL.1Ra > 10, NA, exam1_MESA$IL.1Ra)

vars <- c("IL.1.R.AcP", "IL.6", "IL.1.sRI", "IL.1b", "IL.1Ra")
summary(exam1_MESA[, vars])
for(i in vars){
  hist(exam1_MESA[[i]], main = paste0("Histogram of ", i), xlab = paste0("exam1_MESA$", i))
  boxplot(exam1_MESA[[i]], main = paste0("Boxplot of ", i), xlab = paste0("exam1_MESA$", i))
}

summary(lm(chr2.113083453_G ~ IL.1.R.AcP, data=exam1_MESA))
summary(lm(chr2.113083453_G ~ IL.6, data=exam1_MESA))
summary(lm(chr2.113083453_G ~ IL.1.sRI, data=exam1_MESA))
summary(lm(chr2.113083453_G ~ IL.1b, data=exam1_MESA))
summary(lm(chr2.113083453_G ~ IL.1Ra, data=exam1_MESA))

exam5_MESA <- exam5 %>% inner_join(., IL6_MESA_encore, by = c("IID" = "sample.id")) %>% 
  inner_join(., crp_encore[crp_encore$study.x == "MESA", c("sample.id", "crp")], by = c("IID" = "sample.id")) %>%
  mutate(., IL.1Ra = log(IL.1Ra), IL.1.R.AcP = log(IL.1.R.AcP), IL.6 = log(IL.6), IL.1.sRI = log(IL.1.sRI), IL.1b = log(IL.1b))

# Delete observations with log(IL6) > 8
exam5_MESA$IL.6 <- ifelse(exam5_MESA$IL.6 > 7, NA, exam5_MESA$IL.6)
# Delete observations with log(IL1B) > 10
exam5_MESA$IL.1b <- ifelse(exam5_MESA$IL.1b > 10, NA, exam5_MESA$IL.1b)
# delete log(IL1Ra) > 10
exam5_MESA$IL.1Ra <- ifelse(exam5_MESA$IL.1Ra > 10, NA, exam5_MESA$IL.1Ra)

summary(exam5_MESA[, vars])
for(i in vars){
  hist(exam5_MESA[[i]], main = paste0("Histogram of ", i), xlab = paste0("exam5_MESA$", i))
  boxplot(exam5_MESA[[i]], main = paste0("Boxplot of ", i), xlab = paste0("exam5_MESA$", i))
}

summary(lm(chr2.113083453_G ~ IL.1.R.AcP, data=exam5_MESA))
summary(lm(chr2.113083453_G ~ IL.6, data=exam5_MESA))
summary(lm(chr2.113083453_G ~ IL.1.sRI, data=exam5_MESA))
summary(lm(chr2.113083453_G ~ IL.1b, data=exam5_MESA))
summary(lm(chr2.113083453_G ~ IL.1Ra, data=exam5_MESA))

CHS_lv <- CHS_addtl %>% mutate(., unique_subject_key = paste0("CHS_", Individual_ID)) %>% 
  inner_join(., samples[, c("unique_subject_key", "sample.id")], by = "unique_subject_key") %>%
  inner_join(., rs6734238[, c("IID", "chr2:113083453_G")], by = c("sample.id" = "IID")) %>%
  inner_join(., IL6_CHS_encore, by = "sample.id") %>%
  inner_join(., crp_encore[crp_encore$study.x == "CHS", c("sample.id", "crp")], by = "sample.id") %>%
  mutate(., il1ray2 = log(il1ray2), il1ray5 = log(il1ray5)) %>%
  rename(., IID = sample.id, chr2.113083453_G = `chr2:113083453_G`, IL.1Ra = il1ray2)
CHS_lv2 <- subset(CHS_lv, CHS_lv$il6harmonizedy2 != 0 & CHS_lv$il6harmonizedy2 < mean(CHS_lv$il6harmonizedy2, na.rm = TRUE) + 3*sd(CHS_lv$il6harmonizedy2, na.rm = TRUE)) %>%
  mutate(., il6harmonizedy2 = log(il6harmonizedy2))

## Original IL6 Variable
summary(lm(chr2.113083453_G ~ il6, data=CHS_lv))

## IL6Harmonized2 - after removing outliers and missings
summary(lm(chr2.113083453_G ~ il6harmonizedy2, data=CHS_lv2))

# exam1_analysis <- exam1_MESA[, c("IID", "age", "sex", "race", "bmi1c", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
# exam1_analysis$study <- "MESA_1"
# exam5_analysis <- exam5_MESA[, c("IID", "age", "sex", "race", "bmi5c", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
# exam5_analysis$study <- "MESA_5"
# chs_analysis <- CHS_lv[, c("IID", "age", "sex", "race", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
# chs_analysis$study <- "CHS"
# 
# all_studies <- rbind(exam1_analysis[, -5], exam5_analysis[, -5], chs_analysis)
# 
# chs_analysis_2 <- CHS_lv2[, c("IID", "age", "sex", "race", "chr2.113083453_G", "IL.1Ra", "il6harmonizedy2", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
# chs_analysis_2$study <- "CHS"





