##############################################################################
## Title: CFS Composite Phenotype File Creation for Encore GWAS Analysis
## Version: 1
## Author: Regina Manansala
## Date Created: 31-January-2020
## Date Modified: 25-March-2020
##############################################################################

library(lavaan)
library(data.table)
library(dplyr)
library(stringr)

data_dir <- "~/Documents/Inflammation_Data/"
setwd(data_dir)

## Import Freeze8 Sample Data
samples <- fread(paste0(data_dir,"freeze8_sample_annot_2019-03-28.txt"))

## Get list of sub-directories named for each study included
study <- list.dirs(path = data_dir, full.names = F, recursive = F)

## Import principal component data
pcs <- fread(paste0(data_dir,"pcair_results.txt"), header=F)

## IMPORT INFLAMMATION AND DEMOGRAPHICS DATA ##
for(i in 1:length(study)){
  ## Get files within the study subdirectories
  sub_dir <- list.files(paste0(data_dir, study[i]))
  if(length(grep('topmed_dcc_inflammation_v1', sub_dir)) > 0){
    inf_path <- paste0(data_dir, "/", study[i], "/", sub_dir[grep('topmed_dcc_inflammation_v1', sub_dir)])
    inf_files <- list.files(path = inf_path)
    assign(paste0(study[i], "_inf"), fread(file = paste0(inf_path, "/", inf_files[grep("topmed_dcc_inflammation_v1.txt", inf_files)])))
    dem_path <- paste0(data_dir, "/", study[i], "/", sub_dir[grep('topmed_dcc_demographic_v3', sub_dir)])
    dem_files <- list.files(path = dem_path)
    assign(paste0(study[i], "_dem"), fread(file = paste0(dem_path, "/", dem_files[grep("topmed_dcc_demographic_v3.txt", dem_files)])))
  } 
  if(length(grep('topmed_dcc_inflammation_v1', sub_dir)) == 0){
    if(study[i] != "FHS"){
      assign(paste0(study[i], "_inf"), fread(file = paste0(data_dir, study[i], "/", sub_dir[grep(".csv", sub_dir, ignore.case = T)])))
    }
    if(study[i] == "FHS"){
      assign(paste0(study[i], "_inf"), fread(file = paste0(data_dir, study[i], "/", sub_dir[grep("inflammation.*.csv", sub_dir, ignore.case = T)])))
    }
  }
}

## Update Study List
study <- sub("_inf", "", ls(pattern="_inf"))

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
    # Find duplicate measures and drop based on assay type frequency
    foo2 <- foo2 %>% 
      group_by(CRP_ASSAY) %>% 
      mutate(ASSAY_NUM = table(CRP_ASSAY)) %>%
      ungroup() %>% 
      arrange(dbGaP_ID, ASSAY_NUM) %>%
      subset(CRP_ASSAY != "Latex-enhanced nephelometry (N High Sensitivity CRP assay) on BN II nephelometer (Dade Behring, Inc.)")
    
    # dim(foo2[duplicated(foo2$unique_subject_key) | duplicated(foo2$unique_subject_key, fromLast = TRUE), ])
    foo3 <- foo2[!duplicated(foo2$unique_subject_key, fromLast = TRUE),]
    assign(study[i], foo3)
  }
}
## ADD SD AND RACE EXCLUSION CRITERIA

#WHI 
varlist <- grep("(^crp($|_1)|^il6($|_1)|^il8($|_1)|^il10($|_1)|^il18($|_1)|^icam($|_1)|^tnfa($|_1)|^pselectin($|_1)|^eselectin($|_1)|^l1_beta($|_1)|^tnfa_r1($|_1)|^mmp1($|_1)|^mmp9($|_1)|^cd40($|_1)|^isoprostane_8_epi_pgf2a($|_1)|^lppla2_act($|_1)|^lppla2_mass($|_1)|^mcp1($|_1)|^mpo($|_1)|^opg($|_1)|^tnfr2($|_1))", 
                colnames(WHI), value = TRUE, ignore.case = TRUE)
for(j in 1:length(varlist)){
  # varcol <- get(study[i])[[varlist[j]]] %>% as.numeric()
  if(sum(is.na(WHI[, varlist[j]])) == nrow(WHI)){
    varlist[j] <- NA
  }     
}
varlist <- varlist[!is.na(varlist)]

## Transform inflammation phenotypes
trans.df <- WHI
for(k in 1:length(varlist)){
  if(is.numeric(trans.df[[varlist[k]]]) == FALSE){
    trans.df <- trans.df %>% mutate_at(varlist[k], as.numeric)
  }
  int <- function(x, na.rm = FALSE) (qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  trans.df <- trans.df %>% mutate_at(varlist, int)
}
assign(paste("WHI", "trans", sep = "_"), trans.df)

## Calculate composite phenotype
# mod <- "comp.pheno =~ IL6 + CRP + IL8 + IL10 + ICAM + ESELECTIN + MMP9 + TNFA + TNFA_R1"
# fit <- cfa(mod, data=WHI_trans, missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

# mod <- "comp.pheno =~ IL6 + CRP + IL8 + ICAM + ESELECTIN + MMP9 + TNFA + TNFA_R1"
# fit <- cfa(mod, data=WHI_trans, missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

# mod <- "comp.pheno =~ IL6 + CRP + IL8 + ESELECTIN + MMP9 + TNFA + TNFA_R1"
# fit <- cfa(mod, data=WHI_trans, missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

# mod <- "comp.pheno =~ IL6 + CRP + ESELECTIN + MMP9 + TNFA + TNFA_R1"
# fit <- cfa(mod, data=WHI_trans, missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

mod <- "comp.pheno =~ IL6 + CRP + ESELECTIN + TNFA + TNFA_R1"
fit <- cfa(mod, data=WHI_trans, missing = "ml")
summary(fit, standardized = T, fit.measures = T)

# mod <- "comp.pheno =~ IL6 + CRP + ESELECTIN + TNFA"
# fit <- cfa(mod, data=WHI_trans, missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

# mod <- "comp.pheno =~ IL6 + CRP + ESELECTIN"
# fit <- cfa(mod, data=WHI_trans, missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

fit_id <- fit@Data@case.idx[[1]]
pred <- data.frame(predict(fit), id = fit_id)
# https://groups.google.com/forum/#!msg/lavaan/UPrU8qG5nOs/70OyCU-1u4EJ
WHI_lv <- tibble::rownames_to_column(trans.df, "id") %>% mutate(id = as.numeric(id)) %>% left_join(., pred, by = "id") %>% dplyr::select(-1)

## Format and export
WHI_encore_prelim <- left_join(WHI_lv[, c("unique_subject_key", "comp.pheno", "SEX", "ANCESTRY")], samples[!duplicated(unique_subject_key), c("sample.id", "unique_subject_key")], by = "unique_subject_key")

write.table(WHI_encore_prelim[!is.na(WHI_encore_prelim$comp.pheno),], "../Inflammation_SEM/WHI_encore_prelim.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## Combine WHI and MESA composite phenotypes and export
MESA_lv <- fread("../Inflammation_SEM/Encore_Prelim/MESA_encore_prelim.txt")
names(MESA_lv)[3:4] <- c("SEX", "ANCESTRY")
MESA_WHI <- rbind(MESA_lv, WHI_encore_prelim[!is.na(WHI_encore_prelim$comp.pheno),])

write.table(MESA_WHI, "../Inflammation_SEM/MESA_WHI_encore.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

