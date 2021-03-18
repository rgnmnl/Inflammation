##############################################################################
## Title: IL1RN Initial Correlations
## Version: 1
## Author: Regina Manansala
## Date Created: 26-February-2020
## Date Modified: 07-April-2020
## Results File: Exam1_Analysis.txt; Exam5_Analysis.txt; CHS_Analysis.txt
##############################################################################

library(tidyr)
library(ggplot2)
library(qqman)
library(data.table)
library(dplyr)
library(lavaan)

data_dir <- "~/Documents/Inflammation_SEM/IL6_Analysis/"
setwd(data_dir)

## Load IL1RN files
rs6734238 <- fread("freeze8_rs6734238.raw")
exam1 <- fread("IL1RN_phenotypes_exam1.txt")
exam5 <- fread("IL1RN_phenotypes_exam5.txt")
CHS_addtl <- fread("CHS_MLB_20200121_CHIP_INFLAMMATION.csv")
samples <- fread("~/Documents/Inflammation_Data/freeze8_sample_annot_2019-03-28.txt")

## Correlations between Inflammation Markers and SNP

exam1_pheno <- names(exam1[, c(6:17)])
exam5_pheno <- names(exam5[, c(12:23)])

# test <- lm(chr2.113083453_G ~ IL.1a, data = exam1)

exam1_corr <- matrix(nrow = length(exam1_pheno), ncol = 2)
dimnames(exam1_corr) <- list(exam1_pheno, c("Beta", "P.val"))

for(i in 1:length(exam1_pheno)){
  df <- exam1
  mod <- lm(df$chr2.113083453_G ~ df[[exam1_pheno[i]]])  
  exam1_corr[i, 1] <- summary(mod)$coef[2,1]
  exam1_corr[i, 2] <- summary(mod)$coef[2,4]
}

exam5_corr <- matrix(nrow = length(exam5_pheno), ncol = 2)
dimnames(exam5_corr) <- list(exam5_pheno,  c("Beta", "P.val"))

for(i in 1:length(exam5_pheno)){
  df <- exam5
  mod <- lm(df$chr2.113083453_G ~ df[[exam5_pheno[i]]])  
  exam5_corr[i, 1] <- summary(mod)$coef[2,1]
  exam5_corr[i, 2] <- summary(mod)$coef[2,4]
}

summary(lm(exam5$chr2.113083453_G ~ exam5$IL.6))

chs_pheno <- grep("^il*|^crp*|^ddi*", names(CHS_addtl), value = TRUE)

chs_corr <- matrix(nrow = length(chs_pheno), ncol = 2)
dimnames(chs_corr) <- list(chs_pheno, c("Beta", "P.val"))
df <- CHS_addtl %>% mutate(., unique_subject_key = paste0("CHS_", Individual_ID)) %>% 
  left_join(., samples[, c("unique_subject_key", "sample.id")], by = "unique_subject_key") %>%
  left_join(., rs6734238[, c("IID", "chr2:113083453_G")], by = c("sample.id" = "IID"))
df$rsID <- as.numeric(df$`chr2:113083453_G`)
for(i in 1:length(chs_pheno)){
  mod <- lm(df$rsID ~ df[[chs_pheno[i]]])  
  chs_corr[i, 1] <- summary(mod)$coef[2,1]
  chs_corr[i, 2] <- summary(mod)$coef[2,4]
}

# Merge MESA and CHS IL6 variables from original files and run vs. SNP

exam1_MESA <- fread("~/Documents/Inflammation_Data_Results/IL6/il6_MESA_20190826.txt") %>% inner_join(., exam1, by = c("sample.id" = "IID"))
summary(lm(exam1_MESA$chr2.113083453_G ~ exam1_MESA$il6))
summary(lm(exam1_MESA$il6 ~ exam1_MESA$IL.1Ra + exam1_MESA$chr2.113083453_G))
summary(lm(exam1_MESA$crp1 ~ exam1_MESA$il6 + exam1_MESA$IL.1Ra))

exam5_MESA <- fread("~/Documents/Inflammation_Data_Results/IL6/il6_MESA_20190826.txt") %>% inner_join(., exam5, by = c("sample.id" = "IID"))
summary(lm(exam5_MESA$chr2.113083453_G ~ exam5_MESA$il6))
summary(lm(exam5_MESA$il6 ~ exam5_MESA$IL.1Ra + exam5_MESA$chr2.113083453_G))

df <- CHS_addtl %>% mutate(., unique_subject_key = paste0("CHS_", Individual_ID)) %>% 
  left_join(., samples[, c("unique_subject_key", "sample.id")], by = "unique_subject_key") %>%
  left_join(., rs6734238[, c("IID", "chr2:113083453_G")], by = c("sample.id" = "IID"))
il1RN_CHS <- fread("~/Documents/Inflammation_Data_Results/IL6/il6_CHS_20190826.txt") %>% inner_join(., df[, c("sample.id", "chr2:113083453_G")], by = "sample.id")
summary(lm(il1RN_CHS[["chr2:113083453_G"]] ~ il1RN_CHS$il6))







