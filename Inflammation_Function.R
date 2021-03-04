##############################################################################
## Title: Inflammation Function 
## Purpose: Create phenotype file for use in Encore
## Author: Regina Manansala
## Date Created: 27-September-2019
## Date Modified: 24-February-2020
##############################################################################

library(tidyr)
library(ggplot2)
library(qqman)
library(data.table)
library(dplyr)

encore_file_creation <- function(j){
  
  ## Set home directory containing all data
  ## Studies should be organized in subdirectories to easily locate associated inflammation and demographic data
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
      if(all(is.na(foo[[toupper(j)]])) == "FALSE"){
        foo2 <- foo2 %>% 
          group_by(foo2[[paste0(toupper(j), "_ASSAY")]]) %>% 
          mutate(ASSAY_NUM = table(get(paste0(toupper(j), "_ASSAY")))) %>%
          ungroup()
        if(j %in% c("crp", "CRP")){
          # Remove the following assay type if gathering CRP measures
          foo3 <- foo2[foo2$CRP_ASSAY != "Latex-enhanced nephelometry (N High Sensitivity CRP assay) on BN II nephelometer (Dade Behring, Inc.)", ]
        } else {
          foo3 <- foo2
        }
        foo3 <- foo3[order(foo3$dbGaP_ID, foo3$ASSAY_NUM),]
        foo4 <- foo3[!duplicated(foo3$unique_subject_key, fromLast = TRUE),]
        assign(study[i], foo4)
      } else {
        foo3 <- foo2[order(foo2$dbGaP_ID),]
        foo4 <- foo3[!duplicated(foo3$unique_subject_key, fromLast = TRUE),]
        assign(study[i], foo4)
      }
    }
  }
  
  # Check if desired inflammation variable was included in each study. If measure was not taken for a particular study, remove from further data processing
  for(i in 1:length(study)){
    if(length(grep(j, names(get(study[i])), value = F, ignore.case = T)) == 0){
      study[i] <- NA
    }
  }
  
  ## Subset study list to exclude studies without desired inflammation variable.
  study <- study[!is.na(study)]
  
  ## Keep wanted covariates and rename/reformat columns to match amongst all studies.
  for(i in 1:length(study)){
    col_subset <- get(study[i]) %>% select(grep("(_id|shareid)", colnames(get(study[i])), ignore.case = T, value = T), 
                                           grep("unique_subject_key", colnames(get(study[i])), ignore.case = T, value = T),
                                           grep(paste0("age_at_", j, "|", "AGE_FIRSTVISIT"), colnames(get(study[i])), ignore.case = T, value = T),
                                           grep("^AGE$", colnames(get(study[i])), ignore.case = F, value = T),
                                           grep("sex|male", colnames(get(study[i])), ignore.case = T, value = T),
                                           grep("ancestry|race", colnames(get(study[i])), ignore.case = T, value = T),
                                           grep(paste0("^", j, "$", "|","^", j, "_1$"), colnames(get(study[i])), ignore.case = T, value = T))
    # For tnfa, remove tnfa_r1 columns.
    if(j == "tnfa"){
      col_subset <- select(col_subset, -grep("tnfa_r1", colnames(col_subset), ignore.case = T))
    }
    # Create dummy column if no inflammation measure in study
    if(ncol(col_subset) < 6){
      col_subset$miss_inf <- NA 
    }
    # Rename variable names and recode (by study)
    if(study[i] %in% c("CFS", "COPDGene", "GeneSTAR")){ ##"JHS"
      names(col_subset) <- c("subject_id", "unique_subject_key", "age", "sex", "race", j)
      col_subset <- mutate(col_subset, sex = ifelse(sex == 1, "male", "female"), 
                           race = ifelse(race == 1, "White", "Black"))
    } 
    if(study[i] == "WHI"){
      names(col_subset) <- c("subject_id", "unique_subject_key", "age", "sex", "race", j)
      col_subset <- mutate(col_subset, sex = recode(sex, "Female" = "female"), 
                           race = recode(race, 'African American' = "Black", 
                                         "Asian American" = "Asian", 
                                         "European American" = "White",
                                         "Hispanic American" = "Other",
                                         "Native American" = "AI_AN"))
    }
    if(study[i] %in% c("FHS")){
      names(col_subset) <- c("subject_id", "unique_subject_key", "age", "sex", "race", j)
      col_subset <- mutate(col_subset, sex = ifelse(sex == 1, "male", "female"), 
                           race = recode(race, `1` = "White"))
    }
    if(study[i] %in% c("ARIC", "CARDIA", "CHS", "GENOA", "MESA", "OOA", "SOL")){
      col_subset <- col_subset[,!(names(col_subset) %in% c("SUBJECT_ID.y"))]
      names(col_subset) <- c("subject_id", "unique_subject_key", "age", "sex", "race", j)
    }
    if(study[i] == "JHS"){
      names(col_subset) <- c("subject_id", "unique_subject_key", "age", "sex", "race", j)
      col_subset <- mutate(col_subset, sex = ifelse(sex == 1, "male", "female"), race = recode(race, 'AFR' = "Black"))
    }
    
    # Recoding for all study datasets
    if(nrow(col_subset) > 0){
      ## Recode subject ID variable as character
      if(is.character(col_subset$subject_id) == "FALSE"){
        id_recode <- mutate(col_subset, subject_id = as.character(subject_id))
      } else {
        id_recode <- col_subset
      }
      ## Recode age variable as numeric
      if(is.integer(id_recode$age) == "TRUE"){
        age_recode <- mutate(id_recode, age = as.numeric(age))
      } else {
        age_recode <- id_recode
      }
      ## Inflammation variable recode as numeric
      if(is.character(age_recode[[j]]) == "TRUE"){
        pheno_recode <- age_recode
        pheno_recode[[j]] <- as.numeric(pheno_recode[[j]])
      } else {
        pheno_recode <- age_recode
      }
      pheno_recode$study <- study[i]
      
      # Subset inflammation variable to remove missing values and values that exceed 3 standard deviations of the mean
      sub <- subset(pheno_recode, get(j) != 0 & get(j) < mean(get(j), na.rm = TRUE) + 3*sd(get(j), na.rm = TRUE))
      
      # Create a dataframe labeled with study and inflammation variable
      assign(paste(study[i], j, sep = "_"), sub)
    }
  }
  
  ## Combine all studies into one dataset and merge with principal component data
  inf_all_study <- do.call(rbind, lapply(ls(pattern = paste0("_", j)), function(x) get(x))) %>% 
    left_join(., subset(samples, !duplicated(unique_subject_key)), by="unique_subject_key") %>%
    left_join(., pcs, by=c("sample.id" = "V1"))
  names(inf_all_study)[20:30] <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11")
  
  ## Create a race-study variable
  inf_all_study$race_study <- ifelse(inf_all_study$race == "" | is.na(inf_all_study$race), "", paste(inf_all_study$study.x, inf_all_study$race, sep = "_"))
  
  ## Remove observations with race-study frequencies <50 and observations with any missing values
  inf_encore <- inf_all_study[inf_all_study$race_study %in% names(table(inf_all_study$race_study))[table(inf_all_study$race_study) > 50] & 
                                !is.na(inf_all_study$PC1) &
                                !is.na(inf_all_study$age), c("sample.id", j, "age", "sex", "race_study", "PC1", "PC2", "PC3", 
                                                             "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "study.x")]
  ## Log-transform inflammation variable
  inf_encore[[j]] <- log(inf_encore[[j]])
  return(inf_encore)
}


crp_encore <- encore_file_creation(j = "crp")
cd40_encore <- encore_file_creation(j = "cd40")
icam_encore <- encore_file_creation(j = "icam1")
il1_encore <- encore_file_creation(j = "il1_beta")
il6_encore <- encore_file_creation(j = "il6")
il8_encore <- encore_file_creation(j = "il8")
il10_encore <- encore_file_creation(j = "il10")
il18_encore <- encore_file_creation(j = "il18")
isopro_encore <- encore_file_creation(j = "isoprostane_8_epi_pgf2a")
lppla2_act_encore <- encore_file_creation(j = "lppla2_act")
lppla2_mass_encore <- encore_file_creation(j = "lppla2_mass")
mmp1_encore <- encore_file_creation(j = "mmp1")
mmp9_encore <- encore_file_creation(j = "mmp9")
mcp1_encore <- encore_file_creation(j = "mcp1")
mpo_encore <- encore_file_creation(j = "mpo")
opg_encore <- encore_file_creation(j = "opg")
pselectin_encore <- encore_file_creation(j = "pselectin")
eselectin_encore <- encore_file_creation(j = "eselectin")
tnfa_encore <- encore_file_creation(j = "tnfa")
tnfa_r1_encore <- encore_file_creation(j = "tnfa_r1")
tnfr2_encore <- encore_file_creation(j = "tnfr2")

