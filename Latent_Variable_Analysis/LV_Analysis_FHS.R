##############################################################################
## Title: CFS Composite Phenotype File Creation for Encore GWAS Analysis
## Version: 1
## Author: Regina Manansala
## Date Created: 02-March-2020
## Date Modified: 30-March-2020
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
      assign(paste0(study[i], "_og"), fread(file = paste0(data_dir, study[i], "/", sub_dir[grep("inflammation.*.csv", sub_dir, ignore.case = T)])))
      assign(paste0(study[i], "_add"), fread(file = paste0(data_dir, study[i], "/", sub_dir[grep("additional.*.csv", sub_dir, ignore.case = T)])))
      FHS_inf <- left_join(FHS_og, FHS_add, by = "shareid")
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

#FHS 554 vars
varlist <- grep("(^crp($|_1)|^il6($|_1)|^il8($|_1)|^il10($|_1)|^il18($|_1)|^icam($|_1)|^tnfa($|_1)|^pselectin($|_1)|^eselectin($|_1)|^l1_beta($|_1)|^tnfa_r1($|_1)|^mmp1($|_1)|^mmp9($|_1)|^cd40($|_1)|^isoprostane_8_epi_pgf2a($|_1)|^lppla2_act($|_1)|^lppla2_mass($|_1)|^mcp1($|_1)|^mpo($|_1)|^opg($|_1)|^tnfr2($|_1))", 
                colnames(FHS), value = TRUE, ignore.case = TRUE)
for(j in 1:length(varlist)){
  # varcol <- get(study[i])[[varlist[j]]] %>% as.numeric()
  if(sum(is.na(FHS[, varlist[j]])) == nrow(FHS)){
    varlist[j] <- NA
  }     
}
varlist <- varlist[!is.na(varlist)]

## Transform inflammation phenotypes
trans.df <- FHS
for(k in 1:length(varlist)){
  if(is.numeric(trans.df[[varlist[k]]]) == FALSE){
    trans.df <- trans.df %>% mutate_at(varlist[k], as.numeric)
  }
  int <- function(x, na.rm = FALSE) (qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  trans.df <- trans.df %>% mutate_at(varlist, int) 
}
trans.df <- left_join(trans.df, samples[!duplicated(unique_subject_key), c("sample.id", "unique_subject_key")], by = "unique_subject_key")
assign(paste("FHS", "trans", sep = "_"), trans.df)

# FHS_cc <- FHS_trans[complete.cases(FHS_trans[,c("CRP","IL6","IL18","ICAM","TNFA","PSELECTIN","MMP9")]),]

summary(FHS_trans[, c("IL6", "CRP", "IL18", "ICAM", "TNFA", "PSELECTIN", "MMP9", "OPG", "MCP1", "LPPLA2_ACT", "LPPLA2_MASS")])

## Calculate composite phenotype
# mod <- "comp.pheno =~ CRP + IL6 + IL18 + ICAM + TNFA + PSELECTIN + MMP9 + OPG + MCP1 + LPPLA2_ACT + LPPLA2_MASS"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ CRP + IL6 + IL18 + ICAM + TNFA + PSELECTIN + MMP9 + OPG + MCP1 + LPPLA2_MASS"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ CRP + IL6 + IL18 + ICAM + TNFA + MMP9 + OPG + MCP1 + LPPLA2_MASS"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ CRP + IL6 + IL18 + ICAM + TNFA + MMP9 + OPG + MCP1"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ CRP + IL6 + IL18 + ICAM + TNFA + MMP9 + MCP1"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ CRP + IL6 + IL18 + ICAM + TNFA + MMP9 + OPG"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ IL6 + CRP + IL18 + TNFA + MMP9 + LPPLA2_MASS"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ IL6 + CRP + IL18 + TNFA + MMP9 + LPPLA2_ACT"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

mod <- "comp.pheno =~ IL6 + CRP + IL18 + TNFA + MMP9"
fit <- cfa(mod, data=FHS_trans,missing = "ml")
summary(fit, standardized = T, fit.measures = T)

# mod <- "comp.pheno =~ IL6 + CRP + IL18 + ICAM + PSELECTIN + MMP9"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ IL6 + CRP + IL18 + PSELECTIN + MMP9"
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)
# 
# mod <- "comp.pheno =~ IL6 + CRP + ICAM + PSELECTIN + MMP9"
# # RMSEA                                          0.125
# # Comparative Fit Index (CFI)                    0.847
# fit <- cfa(mod, data=FHS_trans,missing = "ml")
# summary(fit, standardized = T, fit.measures = T)

# fit_id <- fit@Data@case.idx[[1]]
# pred <- data.frame(predict(fit), id = fit_id)
# # https://groups.google.com/forum/#!msg/lavaan/UPrU8qG5nOs/70OyCU-1u4EJ
# FHS_lv <- tibble::rownames_to_column(trans.df, "id") %>% mutate(id = as.numeric(id)) %>% left_join(., pred, by = "id") %>% dplyr::select(-1)

## Format and export
# FHS_encore_prelim <- left_join(FHS_lv[, c("sample.id", "comp.pheno", "SEX", "ANCESTRY")], pcs, by = c("sample.id"="V1"))
# names(FHS_encore_prelim)[5:15] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11")
# write.table(FHS_encore_prelim, "../Inflammation_SEM/FHS_encore_prelim.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

