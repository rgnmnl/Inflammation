##############################################################################
## Title: IL1RN Linear Models
## Version: 2
## Author: Regina Manansala
## Date Created: 07-April-2020
## Date Modified: 15-April-2020
## Results File: lm_il1rn.docx
##				 lm_il6.doc
##				 lm_crp.doc
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
  mutate(., IL.1Ra = log(IL.1Ra))

# crp_plot <- gather(exam1_MESA[, c("crp1", "crp")], "variable", "value", c("crp1", "crp"))
# crp_plot$log_value <- ifelse(crp_plot$variable == "crp", crp_plot$value, log(crp_plot$value))
# ggplot(data=crp_plot, aes(x=log_value, group=variable, fill=variable)) +
#   geom_density(adjust=1.5, alpha=.4) +
#   scale_fill_brewer()

exam5_MESA <- exam5 %>% inner_join(., IL6_MESA_encore, by = c("IID" = "sample.id")) %>% 
  inner_join(., crp_encore[crp_encore$study.x == "MESA", c("sample.id", "crp")], by = c("IID" = "sample.id")) %>%
  mutate(., IL.1Ra = log(IL.1Ra))
# crp_plot <- gather(exam1_MESA[, c("CRP", "crp")], "variable", "value", c("CRP", "crp"))
# crp_plot$log_value <- ifelse(crp_plot$variable == "crp", crp_plot$value, log(crp_plot$value))
# ggplot(data=crp_plot, aes(x=log_value, group=variable, fill=variable)) +
#   geom_density(adjust=1.5, alpha=.4) +
#   scale_fill_brewer()

# samples <- fread("~/Documents/Inflammation_Data/freeze8_sample_annot_2019-03-28.txt")

CHS_lv <- CHS_addtl %>% mutate(., unique_subject_key = paste0("CHS_", Individual_ID)) %>% 
  inner_join(., samples[, c("unique_subject_key", "sample.id")], by = "unique_subject_key") %>%
  inner_join(., rs6734238[, c("IID", "chr2:113083453_G")], by = c("sample.id" = "IID")) %>%
  inner_join(., IL6_CHS_encore, by = "sample.id") %>%
  inner_join(., crp_encore[crp_encore$study.x == "CHS", c("sample.id", "crp")], by = "sample.id") %>%
  mutate(., il1ray2 = log(il1ray2), il1ray5 = log(il1ray5)) %>%
  rename(., IID = sample.id, chr2.113083453_G = `chr2:113083453_G`, IL.1Ra = il1ray2)

# crp_plot <- gather(CHS_lv[, c("crpyr2", "crpyr5", "crpyr9", "crp")], "variable", "value", c("crpyr2", "crpyr5", "crpyr9", "crp"))
# crp_plot$log_value <- ifelse(crp_plot$variable == "crp", crp_plot$value, log(crp_plot$value))
# ggplot(data=crp_plot, aes(x=log_value, group=variable, fill=variable)) +
#   geom_density(adjust=1.5, alpha=.4) +
#   scale_fill_brewer()

exam1_analysis <- exam1_MESA[, c("IID", "age", "sex", "race", "bmi1c", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
exam1_analysis$study <- "MESA_1"
exam5_analysis <- exam5_MESA[, c("IID", "age", "sex", "race", "bmi5c", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
exam5_analysis$study <- "MESA_5"
chs_analysis <- CHS_lv[, c("IID", "age", "sex", "race", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
chs_analysis$study <- "CHS"

all_studies <- rbind(exam1_analysis[, -5], exam5_analysis[, -5], chs_analysis)

###################################################################
############ Regress SNP vs. il6, IL1RA; il6 vs. IL1RA ############ 
###################################################################

summary(lm(chr2.113083453_G ~ il6, data=exam1_analysis))
summary(lm(chr2.113083453_G ~ IL.1Ra, data=exam1_analysis))
summary(lm(IL.1Ra ~ il6, data=exam1_analysis))

summary(lm(chr2.113083453_G ~ il6, data=exam5_analysis))
summary(lm(chr2.113083453_G ~ IL.1Ra, data=exam5_analysis))
summary(lm(IL.1Ra ~ il6, data=exam5_analysis))

summary(lm(chr2.113083453_G ~ il6, data=chs_analysis))
summary(lm(chr2.113083453_G ~ IL.1Ra, data=chs_analysis))
summary(lm(IL.1Ra ~ il6, data=chs_analysis))

###################################################################
####################### Run adjusted models #######################
###################################################################
# 
# ## IL6
# 
# #Model 1: IL6 ~ b0 + b1*SNP + b2*covariates
# 
# mod1 <- lm(il6 ~ chr2.113083453_G + age + sex + race + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = all_studies)
# summary(mod1)
# 
# #Model 2: IL6 ~ a0 + a1*SNP + a2*IL1RN + b3*covariates
# 
# mod2 <- lm(il6 ~ chr2.113083453_G + IL.1Ra + age + sex + race + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = all_studies)
# summary(mod2)
# 
# # Then the direct effect is the estimated coefficient a1:
# summary(mod2)$coef[2]
# 
# # The indirect effetct is estimated as b1 - a1.
# summary(mod1)$coef[2] - summary(mod2)$coef[2]
# 
# ## CRP
#   
# # Model 1: CRP ~ b0 + b1*SNP + b2*covariates
# 
# mod3 <- lm(crp ~ chr2.113083453_G + age + sex + race + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = all_studies)
# summary(mod3)
# 
# # Model 2: CRP ~ a0 + a1*SNP + a2*IL6 + a3*covariates
# 
# mod4 <-  lm(crp ~ chr2.113083453_G + il6 + age + sex + race + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = all_studies)
# summary(mod4)
# 
# # Then the direct effect is the estimated coefficient a1:
# summary(mod4)$coef[2]
# 
# # The indirect effetct is estimated as b1 - a1.
# summary(mod3)$coef[2] - summary(mod4)$coef[2]
# 
# # Estimate the direct and indirect effects the same way.
# # If the indirect effect is large (i.e., much of the effect goes through the mediator (IL6)),
# # then we could further test whether it goes through IL1RN:
# 
# # Model 1: CRP ~ b0 + b1*SNP + b2*covariates
# 
# # mod5 <- lm(crp ~ chr2.113083453_G + age + sex + race + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = all_studies)
# # summary(mod5)
# # 
# # # Model 2: CRP ~ a0 + a1*SNP + a2*IL6 + a3*covariates
# # 
# # mod6 <-  lm(crp ~ chr2.113083453_G + IL.1Ra + age + sex + race + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = all_studies)
# # summary(mod6)
# # 
