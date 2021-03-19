##############################################################################
## Title: LDSC SNP Matching
## Version: 1
## Author: Regina Manansala
## Date Created: 08-July-2020
## Date Modified: 10-December-2020
##############################################################################

for(i in foo){
# 	foo_test <- fread(paste0(i, ".sscore")) %>% subset(., IID %in% pheno_test$IID)
	assign(gsub("\\..*","", i), fread(paste0("/raid-05/SPH/pauer/UKBB_Genotypes/PLINK/LD_SUB_BLACK/",i)))
# 	write.table(foo_test, paste0(i, "_test", ".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
}

flip_snps <- data.table()
for(i in foo){
	foo_test <- fread(paste0("/raid-05/SPH/pauer/UKBB_Genotypes/PLINK/LD_SUB_BLACK/",i)) %>% 
		subset(., V5 == "C" & V6 == "G" | V5 == "G" & V6 == "C" | V5 == "A" & V6 == "T" | V5 == "T" & V6 == "A")	
	flip_snps <- rbind(flip_snps, foo_test)
}
	
	V5 %in% c("A", "T") && V6 %in% c("A", "T")

	assign(gsub("\\..*","", i), fread(paste0("/raid-05/SPH/pauer/UKBB_Genotypes/PLINK/LD_SUB_BLACK/",i)))
# 	write.table(foo_test, paste0(i, "_test", ".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
}

head(cp[cp$Allele1 == "C" & cp$Allele2 == "G" | cp$Allele1 == "G" & cp$Allele2 == "C" | cp$Allele1 == "A" & cp$Allele2 == "T" | cp$Allele1 == "T" & cp$Allele2 == "A",])
head(ra[ra$a1 == "C" & ra$a2 == "G" | ra$a1 == "G" & ra$a2 == "C" | ra$a1 == "A" & ra$a2 == "T" | ra$a1 == "T" & ra$a2 == "A",])head(cp[cp$Allele1 == "C" & cp$Allele2 == "G" | cp$Allele1 == "G" & cp$Allele2 == "C" | cp$Allele1 == "A" & cp$Allele2 == "T" | cp$Allele1 == "T" & cp$Allele2 == "A",])
subset(cp, Allele1 == "C" & Allele2 == "G" | Allele1 == "G" & Allele2 == "C" | Allele1 == "A" & Allele2 == "T" | Allele1 == "T" & Allele2 == "A")

##############################################################################

## IL6 (Using Encore Results)

library(data.table)
library(dplyr)
library(compare)

w_hm3 <- fread("~/Data/manansa2/ldscore/w_hm3.snplist")

comp_pheno <- fread("~/Data/manansa2/lv_encore/comp_pheno_white_results.txt")
topmed_snps <- fread("~/Data/manansa2/lv_encore/TOPMed_SNP_rsID_v2.txt")
cp_snps <- left_join(comp_pheno, topmed_snps[, c("V1", "V2", "V4", "V5", "V9")], by = c("CHR" = "V1", "POS" = "V2", "Allele1" = "V4", "Allele2" = "V5"))
cp_snps <- cp_snps %>% mutate(., SNPID = V9) %>% select(., 1:15)
ref_snps <- fread("~/Data/manansa2/ldscore/w_hm3.snplist")

il6 <- fread("IL6/IL6_WHITE_results.txt")
il6_snps <- left_join(il6, topmed_snps[, c("V1", "V2", "V4", "V5", "V9")], by = c("CHR" = "V1", "POS" = "V2", "Allele1" = "V4", "Allele2" = "V5"))
il6_snps <- il6_snps %>% mutate(., SNPID = V9) %>% select(., 1:14)
cp_il6 <- intersect(il6_snps$SNPID, cp_snps$SNPID)
cp_il6_snps <- intersect(cp_il6, ref_snps$SNP)
cp_il6 <- subset(cp_snps, SNPID %in% cp_il6_snps)
il6_cp <- subset(il6_snps, SNPID %in% cp_il6_snps)
write.table(cp_il6, "cp_il6_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(il6_cp, "il6_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## CRP WBC (Using Encore Results)

wbc <- fread("CRP/wbc_white_results.txt")
wbc_snps <- left_join(wbc, topmed_snp[, c("V1", "V2", "V4", "V5", "V9")], by = c("CHR" = "V1", "POS" = "V2", "Allele1" = "V4", "Allele2" = "V5"))
wbc_snps <- wbc_snps %>% mutate(., SNPID = V9) %>% select(., 1:14)
cp_wbc <- intersect(wbc_snps$SNPID, cp_snps$SNPID)
cp_wbc_snps <- intersect(cp_wbc, ref_snps$SNP)
cp_wbc <- subset(cp_snps, SNPID %in% cp_wbc_snps)
wbc_cp <- subset(wbc_snps, SNPID %in% cp_wbc_snps)
write.table(cp_wbc, "cp_wbc_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(wbc_cp, "wbc_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

wbc <- fread("wbc_ldsc.txt")
cp_wbc <- fread("cp_wbc_ldsc.txt")

cp_wbc_2 <- left_join(cp_wbc, w_hm3, by = c("SNPID"="SNP"))
length(which(cp_wbc_2$Allele1 == cp_wbc_2$A1))
length(which(cp_wbc_2$Allele1 != cp_wbc_2$A1))

length(which(cp_wbc_2$Allele2 == cp_wbc_2$A2))
length(which(cp_wbc_2$Allele2 != cp_wbc_2$A2))

cp_wbc_2$Allele1 <- cp_wbc_2$A1
cp_wbc_2$Allele2 <- cp_wbc_2$A2

write.table(cp_wbc_2[, 1:15], "cp_wbc_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

wbc_2 <- left_join(wbc, w_hm3, by = c("SNPID"="SNP"))
length(which(wbc_2$Allele1 == wbc_2$A1))
length(which(wbc_2$Allele1 != wbc_2$A1))

length(which(wbc_2$Allele2 == wbc_2$A2))
length(which(wbc_2$Allele2 != wbc_2$A2))
wbc_2$Allele2 <- wbc_2$A2
wbc_2$Allele1 <- wbc_2$A1

write.table(wbc_2[, 1:14], "wbc_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

##############################################################################

## Lupus
## Reference GWAS: "bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv"
 
#  lupus <- fread("LUPUS/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv")
#  ldsc_snps <- intersect(cp_snps$SNPID, lupus$rsid) %>% intersect(., ref_snps$SNP)
#  cp_lupus <- subset(cp_snps, SNPID %in% ldsc_snps)
#  lupus_ldsc <- subset(lupus, rsid %in% ldsc_snps)
#  write.table(cp_lupus, "LUPUS/cp_lupus_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(lupus_ldsc, "LUPUS/lupus_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## CAD
## Reference GWAS: "CARDIoGRAM_GWAS_RESULTS.txt" (CAD)
 
#  cad <- fread("CAD/CARDIoGRAM_GWAS_RESULTS.txt")
#  ldsc_snps <- intersect(cp_snps$SNPID, cad$SNP) %>% intersect(., ref_snps$SNP)
#  cp_cad <- subset(cp_snps, SNPID %in% ldsc_snps)
#  cad_ldsc <- subset(cad, SNP %in% ldsc_snps)
#  write.table(cp_cad, "CAD/cp_cad_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(cad_ldsc, "CAD/cad_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  
cad <- fread("cad_ldsc.txt")
cp_cad <- fread("cp_cad_ldsc.txt")
cp_cad_2 <- left_join(cp_cad, w_hm3, by = c("SNPID"="SNP"))
length(which(cp_cad_2$Allele1 == cp_cad_2$A1))
length(which(cp_cad_2$Allele1 != cp_cad_2$A1))
cp_cad_2$Allele1 <- cp_cad_2$A1
cp_cad_2$Allele2 <- cp_cad_2$A2

cad_2 <- left_join(cad, w_hm3, by = "SNP")
length(which(cad_2$reference_allele == cad_2$A1))
length(which(cad_2$reference_allele != cad_2$A1))
cad_2$other_allele <- cad_2$A2
cad_2$reference_allele <- cad_2$A1

write.table(cp_cad_2[, 1:15], "cp_cad_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
write.table(cad_2[, 1:12], "cad_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

##############################################################################

## Crohns Disease
## Reference GWAS: "CD_trans_ethnic_association_summ_stats_b37.txt" (CROHNS DISEASE)   
 
#  crohns <- fread("IBD/CD_trans_ethnic_association_summ_stats_b37.txt")
#  ldsc_snps <- intersect(cp_snps$SNPID, crohns$SNP) %>% intersect(., ref_snps$SNP)
#  cp_crohns <- subset(cp_snps, SNPID %in% ldsc_snps)
#  crohns_ldsc <- subset(crohns, SNP %in% ldsc_snps)
# 
#  crohns_sub <- select(crohns_ldsc, c(1:12, 22)) %>% 
#  	rename(., A1 = A1_effect, A2 = A2_other, N = Mantra_n_samples, BETA = beta_EUR, SE = se_EUR, P = P_EUR, EAF = EAF_EUR)
#  
#  write.table(cp_crohns, "IBD/cp_crohns_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(crohns_sub, "IBD/crohns_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## Schizophrenia                                             
## Reference GWAS: "daner_natgen_pgc_eur" (SCHIZOPHRENIA)    
 
 sch <- fread("SCHIZO/daner_natgen_pgc_eur")
 ldsc_snps <- intersect(cp_snps$SNPID, sch$SNP) %>% intersect(., ref_snps$SNP)
 cp_sch <- subset(cp_snps, SNPID %in% ldsc_snps)
 sch_ldsc <- subset(sch, SNP %in% ldsc_snps)
 write.table(cp_sch, "SCHIZO/cp_sch_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
 write.table(sch_ldsc, "SCHIZO/sch_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
                                                
##############################################################################

## IBD
## Reference GWAS: "IBD_trans_ethnic_association_summ_stats_b37.txt" (INFLAMMATORY BOWEL DISEASE)                                                                                          

#  ibd <- fread("IBD/IBD_trans_ethnic_association_summ_stats_b37.txt")
#  ldsc_snps <- intersect(cp_snps$SNPID, ibd$SNP) %>% intersect(., ref_snps$SNP)
#  cp_ibd <- subset(cp_snps, SNPID %in% ldsc_snps)
#  ibd_ldsc <- subset(ibd, SNP %in% ldsc_snps)
#  
#  ibd_sub <- select(ibd_ldsc, c(1:12, 22)) %>% 
#  	rename(., A1 = A1_effect, A2 = A2_other, N = Mantra_n_samples, BETA = beta_EUR, SE = se_EUR, P = P_EUR, EAF = EAF_EUR)
#  
#  write.table(cp_ibd, "IBD/cp_ibd_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(ibd_sub, "IBD/ibd_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## Stroke
## Reference GWAS: "MEGASTROKE.2.AIS.EUR.out" (STROKE)

#  str <- fread("STROKE/MEGASTROKE.2.AIS.EUR.out")
#  ldsc_snps <- intersect(cp_snps$SNPID, str$MarkerName) %>% intersect(., ref_snps$SNP)
#  cp_str <- subset(cp_snps, SNPID %in% ldsc_snps)
#  str_ldsc <- rename(str, SNP=MarkerName) %>% subset(., SNP %in% ldsc_snps)
#  write.table(cp_str, "STROKE/cp_str_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(str_ldsc, "STROKE/str_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  
#
#                       
# cp_str_2 <- left_join(cp_str, w_hm3, by = c("SNPID" = "SNP"))
# length(which(cp_str_2$Allele1 == cp_str_2$A1))
# length(which(cp_str_2$Allele1 != cp_str_2$A1))
# cp_str_2$Allele1 <- cp_str_2$A1
# cp_str_2$Allele2 <- cp_str_2$A2
#
# 
# str_2 <- left_join(str, w_hm3, by = "SNP")
# length(which(toupper(str_2$Allele1) != str_2$A1))
# length(which(toupper(str_2$Allele1) == str_2$A1))
# str_2$Allele2 <- str_2$A2
# str_2$Allele1 <- str_2$A1
#
# 
# write.table(cp_str_2[, 1:15], "cp_str_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
# write.table(str_2[, 1:7], "str_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

##############################################################################

## Diabetes
## Reference GWAS: "METAANALYSIS_DIAGRAM_SE1.txt" (DIABETES) 
## Note: Need to match chr:pos using liftover and get rsIDs from TopMed - Build 37 to 38    
	
# 	diab <- fread("DIABETES/Diabetes_w_rsID.txt")
# 	ldsc_snps <- intersect(cp_snps$SNPID, diab$SNP) %>% intersect(., ref_snps$SNP)
# 	cp_diab <- subset(cp_snps, SNPID %in% ldsc_snps)
#  	diab_ldsc <- rename(diab, SNP=V9) %>% subset(., SNP %in% ldsc_snps)
#  	write.table(cp_diab, "DIABETES/cp_diab_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  	write.table(diab_ldsc, "DIABETES/diab_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

diab <- fread("diab_ldsc.txt")
cp_diab <- fread("cp_diab_ldsc.txt")

cp_diab_2 <- left_join(cp_diab, w_hm3, by = c("SNPID"="SNP"))
length(which(cp_diab_2$Allele1 == cp_diab_2$A1))
length(which(cp_diab_2$Allele1 != cp_diab_2$A1))
cp_diab_2$Allele1 <- cp_diab_2$A1
cp_diab_2$Allele2 <- cp_diab_2$A2

diab_2 <- left_join(diab, w_hm3, by = "SNP")
length(which(diab_2$Allele1 == diab_2$A1))
length(which(diab_2$Allele1 != diab_2$A1))

length(which(diab_2$Allele2 == diab_2$A2))
length(which(diab_2$Allele2 != diab_2$A2))
diab_2$Allele2 <- diab_2$A2
diab_2$Allele1 <- diab_2$A1

write.table(cp_diab_2[, 1:15], "cp_diab_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
write.table(diab_2[, 1:9], "diab_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

##############################################################################

## UCI
## Reference GWAS: "UC_trans_ethnic_association_summ_stats_b37.txt" (ULCERATIVE COLITIS)    

#  ucol <- fread("IBD/UC_trans_ethnic_association_summ_stats_b37.txt")
#  ldsc_snps <- intersect(cp_snps$SNPID, ucol$SNP) %>% intersect(., ref_snps$SNP)
#  cp_ucol <- subset(cp_snps, SNPID %in% ldsc_snps)
#  ucol_ldsc <- subset(ucol, SNP %in% ldsc_snps)
#  
#  ucol_sub <- select(ucol_ldsc, c(1:12, 22)) %>% 
#  	rename(., A1 = A1_effect, A2 = A2_other, N = Mantra_n_samples, BETA = beta_EUR, SE = se_EUR, P = P_EUR, EAF = EAF_EUR)
#  
#  write.table(cp_ucol, "IBD/cp_ucol_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(ucol_sub, "IBD/ucol_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## Depression
## Reference GWAS: "UKBiobank_broad_12Jan18.txt" (MAJOR DEPRESSIVE DISORDER)

#  dep <- fread("DEPRESSION/UKBiobank_broad_12Jan18.txt")
#  ldsc_snps <- intersect(cp_snps$SNPID, dep$rsid) %>% intersect(., ref_snps$SNP)
#  cp_dep <- subset(cp_snps, SNPID %in% ldsc_snps)
#  dep_ldsc <- subset(dep, rsid %in% ldsc_snps)
#  write.table(cp_dep, "DEPRESSION/cp_dep_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(dep_ldsc, "DEPRESSION/dep_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## WBC
## Reference GWAS: "wbc_N172435_narrow_form.tsv" (WHITE BLOOD CELLS)
 
#  wbc <- fread("WBC/wbc_N172435_narrow_form.tsv")
#  ldsc_snps <- intersect(cp_snps$SNPID, wbc$rsid) %>% intersect(., ref_snps$SNP)
#  cp_wbc <- subset(cp_snps, SNPID %in% ldsc_snps)
#  wbc_ldsc <- subset(wbc, ID_dbSNP49 %in% ldsc_snps) %>%
#  	rename(., SNP = ID_dbSNP49)
#  write.table(cp_wbc, "WBC/cp_wbc_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(wbc_ldsc, "WBC/wbc_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

wbc <- fread("wbc_ldsc.txt")
cp_wbc <- fread("cp_wbc_ldsc.txt")

cp_wbc_2 <- left_join(cp_wbc, w_hm3, by = c("SNPID"="SNP"))
length(which(cp_wbc_2$Allele1 == cp_wbc_2$A1))
length(which(cp_wbc_2$Allele1 != cp_wbc_2$A1))

length(which(cp_wbc_2$Allele2 == cp_wbc_2$A2))
length(which(cp_wbc_2$Allele2 != cp_wbc_2$A2))

cp_wbc_2$Allele1 <- cp_wbc_2$A1
cp_wbc_2$Allele2 <- cp_wbc_2$A2

write.table(cp_wbc_2[, 1:15], "cp_wbc_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

wbc_2 <- left_join(wbc, w_hm3, by = "SNP")
length(which(wbc_2$A1.x == wbc_2$A1.y))
length(which(wbc_2$A1.x != wbc_2$A1.y))

length(which(wbc_2$A2.x == wbc_2$A2.y))
length(which(wbc_2$A2.x != wbc_2$A2.y))
wbc_2$A2.x <- wbc_2$A2.y
wbc_2$A1.x <- wbc_2$A1.y

names(wbc_2)[5:6] <- c("A1", "A2")
write.table(wbc_2[, 1:14], "wbc_ldsc_v2.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

##############################################################################

## Rheumatoid Arthritis
## Reference GWAS: "stahl_etal_2010-RA-EUR.tsv"

#  rhart <- fread("Rheumatoid_Arthritis/stahl_etal_2010-RA-EUR.tsv")
#  ldsc_snps <- intersect(cp_snps$SNPID, rhart$rsid) %>% intersect(., ref_snps$SNP)
#  cp_rhart <- subset(cp_snps, SNPID %in% ldsc_snps)
#  rhart_ldsc <- subset(rhart, rsid %in% ldsc_snps)
#  write.table(cp_rhart, "Rheumatoid_Arthritis/cp_rhart_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(rhart_ldsc, "Rheumatoid_Arthritis/rhart_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## Blood Pressure
## Reference GWAS: "Stahl-et-al_BP_EUR-2018.txt"

#  bp <- fread("BP/Stahl-et-al_BP_EUR-2018.txt")
#  ldsc_snps <- intersect(cp_snps$SNPID, bp$SNP) %>% intersect(., ref_snps$SNP)
#  cp_bp <- subset(cp_snps, SNPID %in% ldsc_snps)
#  bp_ldsc <- subset(bp, SNP %in% ldsc_snps)
#  write.table(cp_bp, "BP/cp_bp_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(bp_ldsc, "BP/bp_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## CRC
## Reference GWAS: "SAIGE_Paper_CRC-GWAS-UKBB.txt.vcf"

#  crc <- fread("CRC/SAIGE_Paper_CRC-GWAS-UKBB.txt.vcf")
#  ldsc_snps <- intersect(cp_snps$SNPID, crc$ID) %>% intersect(., ref_snps$SNP)
#  cp_crc <- subset(cp_snps, SNPID %in% ldsc_snps)
#  crc_ldsc <- subset(crc, ID %in% ldsc_snps)
#  write.table(cp_crc, "CRC/cp_crc_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#  write.table(crc_ldsc, "CRC/crc_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## Breast Cancer
## Reference GWAS: "oncoarray_bcac_public_release_oct17.txt"

 bca <- fread("BC/oncoarray_bcac_public_release_oct17.txt") %>% select(., c(1:6,23:26))
 bca$SNP <- gsub("\\:.*", "", bca$phase3_1kg_id)
 bca_snp <- subset(bca, grepl("rs", phase3_1kg_id))
 bca_snp <- select(bca_snp, c(3:4, 11, 5:10)) %>%
 	rename(., EAF = bcac_gwas_all_eaf_controls, BETA = bcac_gwas_all_beta, SE = bcac_gwas_all_se, P = bcac_gwas_all_P1df)
 
 ldsc_snps <- intersect(cp_snps$SNPID, bca_snp$SNP) %>% intersect(., ref_snps$SNP)
 
 cp_bca <- subset(cp_snps, SNPID %in% ldsc_snps)
 bca_ldsc <- subset(bca_snp, SNP %in% ldsc_snps)
 write.table(cp_bca, "BC/cp_bca_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
 write.table(bca_ldsc, "BC/bca_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## Asthma
## Reference GWAS: "ADULT1_ADULT2_ONSET_ASTHMA.20180716.allchr.assoc.GC"

 ast <- fread("ASTHMA/ADULT1_ADULT2_ONSET_ASTHMA.20180716.allchr.assoc.GC")
 ldsc_snps <- intersect(cp_snps$SNPID, ast$SNP) %>% intersect(., ref_snps$SNP)
 cp_ast <- subset(cp_snps, SNPID %in% ldsc_snps)
 ast_ldsc <- subset(ast, SNP %in% ldsc_snps)
 write.table(cp_ast, "ASTHMA/cp_ast_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
 write.table(ast_ldsc, "ASTHMA/ast_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################################

## Alzheimers
## Reference GWAS: "Kunkle_etal_AD-EUR-2019.txt"

 alz <- fread("ALZ/Kunkle_etal_AD-EUR-2019.txt")
 ldsc_snps <- intersect(cp_snps$SNPID, alz$SNP) %>% intersect(., ref_snps$SNP)
 cp_alz <- subset(cp_snps, SNPID %in% ldsc_snps)
 alz_ldsc <- subset(alz, SNP %in% ldsc_snps)
 write.table(cp_alz, "ALZ/cp_alz_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
 write.table(alz_ldsc, "ALZ/alz_ldsc.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



