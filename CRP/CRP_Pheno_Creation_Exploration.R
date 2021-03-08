##############################################################################
## Title: CRP Phenotype File Creation 
## Purpose: Create CRP phenotype file for use in Encore
## Author: Regina Manansala
## Date Created: 13-June-2019
## Date Modified: 02-July-2019
##############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

## Set home directory containing all data
setwd("~/Documents/Inflammation_Data/")

## Import Freeze8 Sample Data
samples <- fread("freeze8_sample_annot_2019-03-28.txt")

## Get list of sub-directories named for each study included
study <- list.dirs(full.names = F, recursive = F)

## Import principal component data
pcs <- fread("pcair_results.txt", header=F)

## IMPORT INFLAMMATION AND DEMOGRAPHICS DATA ##

for(i in 1:length(study)){
  ## Get files within the study subdirectories
  sub_dir <- list.files(paste0(getwd(), "/", study[i]))
  if(length(grep('topmed_dcc_inflammation_v1', sub_dir)) > 0){
    inf_path <- paste0(getwd(), "/", study[i], "/", sub_dir[grep('topmed_dcc_inflammation_v1', sub_dir)])
    inf_files <- list.files(path = inf_path)
    assign(paste0(study[i], "_inf"), fread(file = paste0(inf_path, "/", inf_files[grep("topmed_dcc_inflammation_v1.txt", inf_files)])))
    dem_path <- paste0(getwd(), "/", study[i], "/", sub_dir[grep('topmed_dcc_demographic_v3', sub_dir)])
    dem_files <- list.files(path = dem_path)
    assign(paste0(study[i], "_dem"), fread(file = paste0(dem_path, "/", dem_files[grep("topmed_dcc_demographic_v3.txt", dem_files)])))
  } else {
    assign(paste0(study[i], "_inf"), fread(file = paste0(study[i], "/", sub_dir[grep(".csv", sub_dir, ignore.case = T)])))
  }
}

## COMBINE INFLAMMATION AND DEMOGRAPHICS DATA ##
## FILTER BY SAMPLE ID ##

inf <- ls(pattern = "_inf")
dem <- ls(pattern = "_dem")
for(i in 1:length(study)){
  ## Remove observations with 'DS' in consent column
  samples <- samples[samples$study %in% c("GeneSTAR", "CFS") | grepl("(-DS|-DS-|DS-)", samples$consent) == FALSE,]
  
  ## Get all studies with separate demographic data. Combine with associated inflammation data using the unique subject key
  if(length(grep(study[i], dem)) > 0){
    inf_dat <- grep(study[i], inf, value = TRUE)
    dem_dat <- grep(study[i], dem, value = TRUE)
    foo <- left_join(get(inf_dat), get(dem_dat), by = "unique_subject_key")
    foo2 <- foo[foo$unique_subject_key %in% samples$unique_subject_key,]
    assign(study[i], foo2)
  } 
  ## Get all studies with inflammation and demographic data combined (except WHI) and create a unique subject key using study name and ID
  if(length(grep(study[i], dem)) == 0 & study[i] != "WHI"){
    foo <- get(paste0(study[i], "_inf"))
    foo$unique_subject_key <- paste(study[i], foo[[grep("(_id|shareid)", colnames(foo), ignore.case = T, value = T)]], sep = "_") #### HOW TO DEAL WITH UCASE/LCASE IN COLNAME
    foo2 <- foo[foo$unique_subject_key %in% samples$unique_subject_key,]
    assign(study[i], foo2)
  }
  
  ## For WHI, create unique subject key using study name and ID
  if(length(grep(study[i], dem)) == 0 & study[i] == "WHI"){
    foo <- get(paste0(study[i], "_inf"))
    foo$unique_subject_key <- paste(study[i], foo[[grep("(_id|shareid)", colnames(foo), ignore.case = T, value = T)]], sep = "_") #### HOW TO DEAL WITH UCASE/LCASE IN COLNAME
    foo2 <- foo[foo$unique_subject_key %in% samples$unique_subject_key,] 
    foo2 <- foo2 %>% 
      group_by(CRP_ASSAY) %>% 
      mutate(CRP_ASSAY_NUM = table(CRP_ASSAY)) %>%
      ungroup()
    foo3 <- foo2[foo2$CRP_ASSAY != "Latex-enhanced nephelometry (N High Sensitivity CRP assay) on BN II nephelometer (Dade Behring, Inc.)", ]
    foo3 <- foo3[order(foo3$dbGaP_ID, foo3$CRP_ASSAY_NUM),]
    foo4 <- foo3[!duplicated(foo3$unique_subject_key, fromLast = TRUE),]
    assign(study[i], foo4)
  }
}

## Keep wanted covariates and rename/reformat columns to match amongst all studies.
for(i in 1:length(study)){
  crp <- get(study[i]) %>% select(grep("(_id|shareid)", colnames(get(study[i])), ignore.case = T, value = T), 
                              grep("unique_subject_key", colnames(get(study[i])), ignore.case = T, value = T),
                              grep("age_at_crp|AGE_FIRSTVISIT", colnames(get(study[i])), ignore.case = T, value = T),
                              grep("^AGE$", colnames(get(study[i])), ignore.case = F, value = T),
                              grep("sex|male", colnames(get(study[i])), ignore.case = T, value = T),
                              grep("ancestry|race", colnames(get(study[i])), ignore.case = T, value = T),
                              grep("crp", colnames(get(study[i])), ignore.case = T, value = T))
  #crp <- setdiff(names(crp), c("SUBJECT_ID.y", "CRP_UNITS", "CRP_ASSAY", "BMI_AT_CRP","SMOKING_AT_CRP"))
  if(study[i] %in% c("CFS", "COPDGene", "GeneSTAR")){ ##"JHS"
    crp <- select(crp, -c("CRP_UNITS", "CRP_ASSAY", "BMI_AT_CRP","SMOKING_AT_CRP"))
  } 
  if(study[i] == "WHI"){
    crp <- select(crp, -c("CRP_UNITS", "CRP_ASSAY", "BMI_AT_CRP","SMOKING_AT_CRP", "CRP_ASSAY_NUM"))
  }
  if(study[i] == "FHS"){
    crp <- select(crp, -c("CRP_UNITS", "CRP_ASSAY"))
  }
  if(study[i] %in% c("ARIC", "CARDIA", "CHS", "GENOA", "MESA", "OOA", "SOL")){
    crp <- crp[,!(names(crp) %in% c("SUBJECT_ID.y"))]
  }
  if(nrow(crp) > 0){
    assign(paste0(study[i], "_crp"), crp)
  }
}

rm(list=ls(pattern = "_inf"))
rm(list=ls(pattern = "_dem"))
rm(list=ls(pattern="foo"))

## Remove Genoa from study list
study <- study[-grep("GENOA", study, value = F)]
#View(select(JHS, grep("AGE", names(JHS), ignore.case = T, value = T)))

## Combine all STUDY_crp data sets into list and iterate over each df in list and rename all variables
## Unlist dfs
#crp_names <- do.call(rbind, lapply(mget(ls(pattern = "_crp")), function(x) colnames(x)))
crp_list <- setNames(lapply(ls(pattern = "_crp"), function(x) get(x)), ls(pattern = "_crp"))
crp_rename <- lapply(crp_list, setNames, c("subject_id", "unique_subject_key", "age", "sex", "race", "crp"))
list2env(crp_rename, .GlobalEnv)

## Get summary of CRP values for each study
for(i in 1:length(study)){
  crp <- get(paste0(study[i], "_crp"))$crp
  print(study[i])
  print(summary(crp))
}

## ID Recoding
for(i in 1:length(study)){
  if(is.character(get(paste0(study[i], "_crp"))$subject_id) == "FALSE"){
    id_recode <- mutate(get(paste0(study[i], "_crp")), subject_id = as.character(subject_id))
  } else {
    id_recode <- get(paste0(study[i], "_crp"))
  }
  assign(paste0(study[i], "_crp"), id_recode)
}

## AGE Recoding
for(i in 1:length(study)){
  if(is.integer(get(paste0(study[i], "_crp"))$age) == "TRUE"){
    age_recode <- mutate(get(paste0(study[i], "_crp")), age = as.numeric(age))
  } else {
    age_recode <- get(paste0(study[i], "_crp"))
  }
  age_recode$study <- study[i]
  assign(paste0(study[i], "_crp"), age_recode)
}

## SEX Recoding
CFS_crp <- mutate(CFS_crp, sex = ifelse(sex == 1, "male", "female"))
COPDGene_crp <- mutate(COPDGene_crp, sex = ifelse(sex == 1, "male", "female"))
FHS_crp <- mutate(FHS_crp, sex = ifelse(sex == 1, "male", "female"))
GeneSTAR_crp <- mutate(GeneSTAR_crp, sex = ifelse(sex == 1, "male", "female"))
JHS_crp <- mutate(JHS_crp, sex = ifelse(sex == 1, "male", "female"))
WHI_crp <- mutate(WHI_crp, sex = recode(sex, "Female" = "female"))

## RACE Recoding
#1 = European American	2 = African American	3 = Asian American	4 = Hispanic American	5 = Native American	6 = Other
CFS_crp <- mutate(CFS_crp, race = ifelse(race == 1, "White", "Black"))
#ifelse(CFS_crp$race == "1",, CFS_crp$race == "Black")
COPDGene_crp <- mutate(COPDGene_crp, race = ifelse(race == 1, "White", "Black"))
FHS_crp <- mutate(FHS_crp, race = recode(race, `1` = "White"))
GeneSTAR_crp <- mutate(GeneSTAR_crp, race = ifelse(race == 1, "White", "Black"))
JHS_crp <- mutate(JHS_crp, race = recode(race, 'AFR' = "Black"))
WHI_crp <- mutate(WHI_crp, race = recode(race, 'African American' = "Black", "Asian American" = "Asian", 
                                         "European American" = "White", "Hispanic American" = "Other",
                                         "Native American" = "AI_AN"))

## CRP Recoding
CFS_crp$crp <- as.numeric(CFS_crp$crp)
for(i in 1:length(study)){
  # if(study[i] == "CFS"){
  #   crp_sub <- mutate(get(paste0(study[i], "_crp")), crp = as.numeric(crp))
  #   crp_sub <- subset(get(paste0(study[i], "_crp")), crp != 0 & crp < mean(crp, na.rm = TRUE) + 3*sd(crp, na.rm = TRUE))
  #   assign(paste0(study[i], "_crp_sub"), crp_sub)
  # } else {
    crp_sub <- subset(get(paste0(study[i], "_crp")), crp != 0 & crp < mean(crp, na.rm = TRUE) + 3*sd(crp, na.rm = TRUE))
    assign(paste0(study[i], "_crp_sub"), crp_sub)
  # }
}
#GeneSTAR_crp$crp_mgL <- GeneSTAR_crp$crp / 1000

## Combine Data
crp <- do.call(rbind, lapply(ls(pattern = "_crp_sub"), function(x) get(x))) %>% 
  left_join(., subset(samples, !duplicated(unique_subject_key)), by="unique_subject_key") %>%
  left_join(., pcs, by=c("sample.id" = "V1"))
names(crp)[20:30] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11")

crp$race_study <- ifelse(crp$race == "" | is.na(crp$race), "", paste(crp$study.x, crp$race, sep = "_"))

crp_encore <- crp[!is.na(crp$PC1) & crp$study.x != "GeneSTAR" & 
                  !(crp$race_study %in% c("", "CHS_AI_AN", "CHS_Asian", "CHS_Other", "COPDGene_Black", "WHI_AI_AN", "WHI_Asian")) & 
                  !is.na(crp$age), 
                  c("sample.id", "crp", "age", "sex", "race_study", "PC1", "PC2", "PC3", 
                    "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
crp_encore$crp <- log(crp_encore$crp)
# write.table(crp_encore, "log_crp_encore.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# write.table(crp_encore[,-5], "norace_crp_encore.txt", sep = "\t", row.names = F, col.names = T, quote = F)

## Run test model for plotting
foo <- lm(crp ~ age + sex + race_study + PC1 + PC2 + PC3 + PC4 + 
  PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=crp_encore)
summary(foo)

## Boxplots for Each Study
explore_plots <- function(df, pheno_col, pheno_name){
  boxplot_by_study <- ggplot(df, aes(x=study.x, y=pheno_col)) + #, color=aims_subtype)) +
    geom_boxplot(outlier.size=.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=5, angle = 60, hjust = 1)) +
    ggtitle(paste0(pheno_name, " Distributions")) +
    xlab("Study") +
    ylab(pheno_name)
  
  log_boxplot_by_study <- ggplot(df, aes(x=study.x, y=log(pheno_col))) +
    geom_boxplot(outlier.size=.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=5, angle = 60, hjust = 1)) +
    ggtitle(paste0("log ", pheno_name, " Distributions")) +
    xlab("Study") +
    ylab(paste0("log ", pheno_name))
  
  density_by_study <- ggplot(df, aes(pheno_col, group=study.x)) +
    facet_wrap(~ study.x, scales = "free") +
    geom_density() +
    theme_minimal() +
    ggtitle(study[i]) +
    xlab(pheno_name)
    
  log_density_by_study <- ggplot(df, aes(log(pheno_col), group=study.x)) +
    facet_wrap(~ study.x, scales = "free") +
    geom_density() +
    theme_minimal() +
    ggtitle(study[i]) +
    xlab(paste0("log ", pheno_name))
  
  return(list(boxplot_by_study = boxplot_by_study, log_boxplot_by_study = log_boxplot_by_study,
              density_by_study = density_by_study, log_density_by_study = log_density_by_study))
}

print(explore_plots(crp, crp$crp, "CRP")$boxplot_by_study)
print(explore_plots(crp, crp$crp, "CRP")$log_boxplot_by_study)
print(explore_plots(crp, crp$crp, "CRP")$density_by_study)
print(explore_plots(crp, crp$crp, "CRP")$log_density_by_study)

## EXPLORATORY PLOTS
# ggplot(crp[crp$study.x == "GeneSTAR", ], aes(crp, color = topmed_project, group=topmed_project)) +
#   #facet_wrap(~ study.x, scales = "free") +
#   geom_density() +
#   theme_minimal() 
# 
# ggplot(GeneSTAR, aes(CRP, color = CRP_ASSAY, group=CRP_ASSAY)) +
#   #facet_wrap(~ study.x, scales = "free") +
#   geom_density() +
#   theme_minimal() 
# ggplot(GeneSTAR, aes(log(CRP), color = CRP_ASSAY, group=CRP_ASSAY)) +
#   #facet_wrap(~ study.x, scales = "free") +
#   geom_density() +
#   theme_minimal() 
# # +
# #   ggtitle("GeneSTAR CRP by ASSAY") +
# #   xlab("log CRP")
# 
# ggplot(GeneSTAR_crp, aes(y=log(crp/1000))) + #, color=aims_subtype)) +
#   geom_boxplot(outlier.size=.5) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size=5, angle = 60, hjust = 1)) +
#   ggtitle(paste0(pheno_name, " Distributions")) +
#   xlab("Study") +
#   ylab(pheno_name)

##PC PLOTS
pca <- crp[,c("race", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11")]
pca$Color[pca$race == "White"] <- "blue"
pca$Color[pca$race == "Black"] <- "black"
pca$Color[pca$race == "AI_AN"] <- "orange"
pca$Color[pca$race == "Asian"] <- "red"
pca$Color[pca$race == "Other"] <- "green"

par(mfrow = c(3,3))
for(i in 2:10){
  #jpeg(paste0(i, j,".jpg",sep=""), height=5, width=5, units='in', quality=100, res=500)
  plot(pca[,c(i,i+1)], pch=20, cex=0.4, col=pca$Color, main = "Principal Components")
  #dev.off()
}

## Models
par(mfrow = c(4,3))
for(i in 1:length(study)){
  crp_sub <- crp[crp$study.x == study[i] & !is.na(crp$crp),]
  if(length(table(crp_sub$race)) > 1 & study[i] != "WHI"){
    mod <- lm(log(crp) ~ age + sex + race + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=crp_sub)
  } 
  if(length(table(crp_sub$race)) == 1 & study[i] != "JHS"){
    mod <- lm(log(crp) ~ age + sex + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=crp_sub)
  }
  if(study[i] == "JHS"){
    crp_sub <- mutate(crp_sub, crp = ifelse(crp == 0, crp+1, crp))
    mod <- lm(log(crp) ~ age + sex + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=crp_sub)
  }
  if(study[i] == "WHI"){
    mod <- lm(log(crp) ~ age + race + PC1 + PC2 + PC3 + PC4 + 
                PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=crp_sub)
  }
  plot(density(resid(mod)), main = paste0(study[i], " Residual Plot"))
}

# whi_assay <- list(unique(WHI$CRP_ASSAY)) %>% unlist()
# whi_assay <- whi_assay[-3]
# datalist <- list()
# for(i in 1:length(whi_assay)){
#   assay_name <- whi_assay[i]
#   #  <- mean(WHI[WHI$CRP_ASSAY == whi_assay[i], "CRP"], na.rm=TRUE)
#   assay_var <- var(WHI[WHI$CRP_ASSAY == whi_assay[i], "CRP"], na.rm=TRUE)
#   assay_mean <- summary(WHI[WHI$CRP_ASSAY == whi_assay[i], "CRP"])[4]
#   sum_vec <- c(assay_name, assay_mean, assay_var)
#   datalist[[i]] <- sum_vec
# }
# assay_summary = do.call(rbind, datalist)
# write.table(assay_summary, "WHI_assay_summary.txt", sep = "\t", row.names = F, col.names = F, quote = F)


mod_crp_jhs <- lm(log(crp) ~ age + sex + PC1 + PC2 + PC3 + PC4 + 
                     PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=crp[crp$study.x == "JHS",])
mod_crp_ooa.res  <- resid(mod_crp_ooa)
plot(density(mod_crp_ooa.res))
#abline(0, 0)
