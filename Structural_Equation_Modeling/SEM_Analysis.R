##############################################################################
## Title: IL1RN SEM
## Note: Results from this were not used. SEM was ultimately conducted by
##		Helen using MPLUS
## Version: 1
## Author: Regina Manansala
## Date Created: 07-April-2020
## Date Modified: 08-April-2020
## Results File: SEM_Output.docx
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

## Make data for use with lavaan SEM
exam1_analysis <- exam1_MESA[, c("IID", "age", "sex", "race", "bmi1c", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
exam1_analysis$study <- "MESA_1"
exam5_analysis <- exam5_MESA[, c("IID", "age", "sex", "race", "bmi5c", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
exam5_analysis$study <- "MESA_5"
chs_analysis <- CHS_lv[, c("IID", "age", "sex", "race", "chr2.113083453_G", "IL.1Ra", "il6", "crp", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
chs_analysis$study <- "CHS"
all_studies <- rbind(exam1_analysis[, -5], exam5_analysis[, -5], chs_analysis)

# write.table(exam1_analysis, "../Inflammation_SEM/Exam1_Analysis.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(exam5_analysis, "../Inflammation_SEM/Exam5_Analysis.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# write.table(chs_analysis, "../Inflammation_SEM/CHS_Analysis.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

model <- 'crp ~ a*il6 + b*IL.1Ra + c*chr2.113083453_G
il6 ~ d*IL.1Ra + f*chr2.113083453_G
IL.1Ra ~ e*chr2.113083453_G
varthruil6il1ra := a*d*e
varthruil6 := a*f
varthruil1ra := b*e
totind := varthruil6il1ra + varthruil6 + varthruil1ra
totthruil6 := varthruil6il1ra + varthruil6
varonil6 := e*d
totonil6 := varonil6 + f
'
model.fit <- sem(model, data = exam1_analysis, missing = 'fiml.x')
summary(model.fit,  standardized = T, fit.measures = T)

model <- 'crp ~ a*il6 + b*IL.1Ra + c*chr2.113083453_G
il6 ~ d*IL.1Ra + f*chr2.113083453_G
IL.1Ra ~ e*chr2.113083453_G
varthruil6il1ra := a*d*e
varthruil6 := a*f
varthruil1ra := b*e
totind := varthruil6il1ra + varthruil6 + varthruil1ra
totthruil6 := varthruil6il1ra + varthruil6
varonil6 := e*d
totonil6 := varonil6 + f
'
model.fit <- sem(model, data = exam5_analysis, missing = 'fiml.x')
summary(model.fit,  standardized = T, fit.measures = T)

model <- 'crp ~ a*il6 + b*IL.1Ra + c*chr2.113083453_G
il6 ~ d*IL.1Ra + f*chr2.113083453_G
IL.1Ra ~ e*chr2.113083453_G
varthruil6il1ra := a*d*e
varthruil6 := a*f
varthruil1ra := b*e
totind := varthruil6il1ra + varthruil6 + varthruil1ra
totthruil6 := varthruil6il1ra + varthruil6
varonil6 := e*d
totonil6 := varonil6 + f
'
model.fit <- sem(model, data = chs_analysis, missing = 'fiml.x')
summary(model.fit,  standardized = T, fit.measures = T)

model <- 'crp ~ a*il6 + b*IL.1Ra + c*chr2.113083453_G
crp ~~ il6
il6 ~ d*IL.1Ra + f*chr2.113083453_G
IL.1Ra ~ e*chr2.113083453_G
varthruil6il1ra := a*d*e
varthruil6 := a*f
varthruil1ra := b*e
totind := varthruil6il1ra + varthruil6 + varthruil1ra
totthruil6 := varthruil6il1ra + varthruil6
varonil6 := e*d
totonil6 := varonil6 + f
'
model.fit <- sem(model, data = all_studies, missing = 'fiml.x')
summary(model.fit,  standardized = T, fit.measures = T)

