##############################################################################
## Title: LDScore Command Line Codes
## Version: 1
## Author: Regina Manansala
## Date Created: 17-June-2020
## Date Modified: 10-December-2020
##############################################################################

source ~/Data/conda_install/bin/activate ldsc
# ldsc/ldsc.py -h
# ldsc/munge_sumstats.py -h

## Rheumatoid Arthritis

~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ../lv_encore/Rheumatoid_Arthritis/cp_snps_ldsc.txt \
	--N 1178609 \
	--merge-alleles w_hm3.snplist \
	--out comp_pheno_merge

~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ../lv_encore/Rheumatoid_Arthritis/ra_snps_ldsc.txt \
	--N 1178524 \
	--merge-alleles w_hm3.snplist \
	--out rh_art_merge

# LD Score Regression
~/Data/manansa2/ldscore/ldsc/ldsc.py \
	--rg ~/Data/manansa2/ldscore/rh_art_merge.sumstats.gz,~/Data/manansa2/ldscore/comp_pheno_merge.sumstats.gz \
	--ref-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--w-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--samp-prev 0.25,nan \
	--pop-prev 0.01,nan \
	--out ~/Data/manansa2/lv_encore/Rheumatoid_Arthritis/ra_cp
less ra_cp.log

# ldsc/ldsc.py --h2 rh_art.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ra_bip
# less ra_bip.log
# 
# ldsc/ldsc.py --h2 comp_pheno.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out cp_bip
# less cp_bip.log

##############################################################################

## CRP and IL6 Munge Stats

source ~/Data/conda_install/bin/activate ldsc

~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/CRP/cp_crp_ldsc_v2.txt \
	--N 1207289 \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/CRP/cp_crp_merge_v2
	
~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/CRP/crp_ldsc_v2.txt \
	--N 1207285 \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/CRP/crp_merge_v2

~/Data/manansa2/ldscore/ldsc/ldsc.py \
	--rg cp_crp_merge_v2.sumstats.gz,crp_merge_v2.sumstats.gz \
	--ref-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--w-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--out crp_comp_pheno_ldsc_v2
less crp_comp_pheno_ldsc.log

############

~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/IL6/cp_il6_ldsc.txt \
	--N 1207523 \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/IL6/cp_il6_merge
	
~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/IL6/il6_ldsc.txt \
	--N 1207519 \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/IL6/il6_merge

~/Data/manansa2/ldscore/ldsc/ldsc.py \
	--rg cp_il6_merge.sumstats.gz,il6_merge.sumstats.gz \
	--ref-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--w-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--out ~/Data/manansa2/lv_encore/IL6/il6_comp_pheno_ldsc

##############################################################################

## LUPUS

#!/usr/bin/env bash

### 

#SBATCH --job-name=ldsc_munge
#SBATCH --mem-per-cpu=60000

source ~/Data/conda_install/bin/activate ldsc

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4668589/
### Our study comprised 7,219 cases and 15,991 controls of European ancestry

~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/LUPUS/cp_lupus_ldsc.txt \
	--N 1175280 \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/LUPUS/cp_lupus_merge
	
~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/LUPUS/lupus_ldsc.txt \
	--N 1175200 \
	--ignore OR,OR_lower,OR_upper \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/LUPUS/lupus_merge

~/Data/manansa2/ldscore/ldsc/ldsc.py \
	--rg cp_lupus_merge.sumstats.gz,lupus_merge.sumstats.gz \
	--ref-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--w-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--samp-prev nan,0.3 \
	--pop-prev nan,0.01 \
	--out ~/Data/manansa2/lv_encore/LUPUS/lupus_comp_pheno_ldsc

less lupus_comp_pheno_ldsc.log

conda deactivate

##############################################################################

## IBD

## CROHNS

 ## https://www.ibdgenetics.org/projects.html
 ### a total of 6,333 patients with Crohn's disease; 15,056 healthy controls
 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/IBD/cp_crohns_ldsc.txt \
# 	--N 53274 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/IBD/cp_crohns_merge
# 	
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/IBD/crohns_ldsc.txt \
# 	--N 53271 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/IBD/crohns_merge
# 
# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_crohns_merge.sumstats.gz,crohns_merge.sumstats.gz \
# 	--ref-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.3 \
# 	--pop-prev nan,0.01 \
# 	--out ~/Data/manansa2/lv_encore/IBD/crohns_comp_pheno_ldsc
# 
# less crohns_comp_pheno_ldsc.log

##############################################################################

## IBD

 ## https://www.ibdgenetics.org/projects.html
 ### 15,000 patients with Crohn's disease, 12,000 with UC and 21,000 healthy controls

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/IBD/cp_ibd_ldsc.txt \
# 	--N 53290 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/IBD/cp_ibd_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/IBD/ibd_ldsc.txt \
# 	--N 53287 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/IBD/ibd_merge
# 
# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_ibd_merge.sumstats.gz,ibd_merge.sumstats.gz \
# 	--ref-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.56 \
# 	--pop-prev nan,0.01 \
# 	--out ~/Data/manansa2/lv_encore/IBD/ibd_comp_pheno_ldsc

##############################################################################

## ULCERATIVE COLITIS

 ## https://www.ibdgenetics.org/projects.html
 ### a total of 6,782 patients with UC and 20,099 healthy controls
 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/IBD/cp_ucol_ldsc.txt \
# 	--N 53284 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/IBD/cp_ucol_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/IBD/ucol_ldsc.txt \
# 	--N 53281 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/IBD/ucol_merge
# 
# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_ucol_merge.sumstats.gz,ucol_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.25 \
# 	--pop-prev nan,0.01 \
# 	--out ~/Data/manansa2/lv_encore/IBD/ucol_comp_pheno_ldsc

##############################################################################

## CAD

## http://www.cardiogramplusc4d.org/data-downloads/
### CARDIoGRAM GWAS is a meta-analysis of 22 GWAS studies of European descent imputed to HapMap 2 involving 
### 22,233 cases and 64,762 controls - data as published in: 
### Schunkert H, König IR, Kathiresan S, Reilly MP, Assimes TL, Holm H et al. Large-scale association analysis identifies 13 new susceptibility loci for coronary artery disease. 
### Nat Genet. 2011 43: 333-338

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/CAD/cp_cad_ldsc.txt \
# 	--N 1054123 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/CAD/cp_cad_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/CAD/cad_ldsc.txt \
# 	--N 1054064 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/CAD/cad_merge
# 
# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_cad_merge.sumstats.gz,cad_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.26 \
# 	--pop-prev nan,0.06 \
# 	--out ~/Data/manansa2/lv_encore/CAD/cad_comp_pheno_ldsc
	
##############################################################################

## DEPRESSION

## https://www.nature.com/articles/s41467-018-03819-3
### Broad depression (cases = 113,769, controls = 208,811)
	
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/DEPRESSION/cp_dep_ldsc.txt \
# 	--N 1169837 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/DEPRESSION/cp_dep_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/DEPRESSION/dep_ldsc.txt \
# 	--N 1170412 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/DEPRESSION/dep_merge
# 
# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_dep_merge.sumstats.gz,dep_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.35 \
# 	--pop-prev nan,0.1 \
# 	--out ~/Data/manansa2/lv_encore/DEPRESSION/dep_comp_pheno_ldsc

##############################################################################

## SCHIZO

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/
### 11,260 cases and 24,542 controls

~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/SCHIZO/cp_sch_ldsc.txt \
	--N 1175648 \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/SCHIZO/cp_sch_merge

~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
	--sumstats ~/Data/manansa2/lv_encore/SCHIZO/sch_ldsc.txt \
	--N 1175568 \
	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
	--out ~/Data/manansa2/lv_encore/SCHIZO/sch_merge

~/Data/manansa2/ldscore/ldsc/ldsc.py \
	--rg cp_sch_merge.sumstats.gz,sch_merge.sumstats.gz \
	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--samp-prev nan,0.31 \
	--pop-prev nan,0.004 \
	--out ~/Data/manansa2/lv_encore/SCHIZO/sch_comp_pheno_ldsc

##############################################################################

## STROKE

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/STROKE/cp_str_ldsc.txt \
# 	--N 1184761 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/STROKE/cp_str_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/STROKE/str_ldsc.txt \
# 	--N 1184677 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/STROKE/str_merge
# 
~/Data/manansa2/ldscore/ldsc/ldsc.py \
	--rg cp_str_merge_v2.sumstats.gz,str_merge_v2.sumstats.gz \
	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--samp-prev nan,0.31 \
	--pop-prev nan,0.004 \
	--out ~/Data/manansa2/lv_encore/STROKE/str_comp_pheno_ldsc_v2

##############################################################################

## DIABETES

## https://diabetes.diabetesjournals.org/content/66/11/2888
### 26,676 T2D case and 132,532 control subjects of European ancestry
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/DIABETES/cp_diab_ldsc.txt \
# 	--N 1203190 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/DIABETES/cp_diab_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/DIABETES/diab_ldsc.txt \
# 	--N 1203087 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/DIABETES/diab_merge
# 
# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_diab_merge.sumstats.gz,diab_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.17 \
# 	--pop-prev nan,0.01 \
# 	--out ~/Data/manansa2/lv_encore/DIABETES/diab_comp_pheno_ldsc

##############################################################################

## WBC

## https://pubmed.ncbi.nlm.nih.gov/27863252/

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/WBC/cp_wbc_ldsc.txt \
# 	--N 1204540 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/WBC/cp_wbc_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/WBC/wbc_ldsc.txt \
# 	--N 1204515 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/WBC/wbc_merge
# 
~/Data/manansa2/ldscore/ldsc/ldsc.py \
	--rg cp_wbc_merge_v2.sumstats.gz,wbc_merge_v2.sumstats.gz \
	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
	--out ~/Data/manansa2/lv_encore/WBC/wbc_comp_pheno_ldsc_v2
	# 	--samp-prev nan,### \
# 	--pop-prev nan,### \

##############################################################################
	
## Rheumatoid Arthritis #2

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4243840/
### 5,539 autoantibody positive RA cases and 20,169 controls of European descent

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/Rheumatoid_Arthritis/cp_rhart_ldsc.txt \
# 	--N 1115411 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out cp_stahl_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/Rheumatoid_Arthritis/rhart_ldsc.txt \
# 	--N 1115348 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out rhart_stahl_merge

# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg rhart_stahl_merge.sumstats.gz,cp_stahl_merge.sumstats.gz \
# 	--ref-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev 0.22,nan \
# 	--pop-prev 0.01,nan \
# 	--out ~/Data/manansa2/lv_encore/Rheumatoid_Arthritis/rhart_cp_stahl

##############################################################################

## BP

## https://pubmed.ncbi.nlm.nih.gov/31043756/
### 20,352 cases and 31,358 controls of European descent

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/BP/cp_bp_ldsc.txt \
# 	--N 1205623 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/BP/cp_bp_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/BP/bp_ldsc.txt \
# 	--N 1205520 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/BP/bp_merge

# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_bp_merge.sumstats.gz,bp_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.39 \
# 	--pop-prev nan,0.028 \
# 	--out ~/Data/manansa2/lv_encore/BP/bp_comp_pheno_ldsc
	
##############################################################################

## CRC

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6517433/
### 34,627 CRC cases and 71,379 controls of European ancestry

#!/usr/bin/env bash

### 

#SBATCH --job-name=ldsc_munge
#SBATCH --mem-per-cpu=60000

source ~/Data/conda_install/bin/activate ldsc

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/CRC/cp_crc_ldsc.txt \
# 	--N 1207439 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/CRC/cp_crc_merge

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/CRC/crc_ldsc.txt \
# 	--N 1212142 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/CRC/crc_merge

# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_crc_merge.sumstats.gz,crc_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.33 \
# 	--pop-prev nan,0.04 \
	--out ~/Data/manansa2/lv_encore/CRC/crc_comp_pheno_ldsc

##############################################################################

## BC

## https://pubmed.ncbi.nlm.nih.gov/29059683/
### 122,977 cases and 105,974 controls of European ancestry

#!/usr/bin/env bash

### 

#SBATCH --job-name=ldsc_munge
#SBATCH --mem-per-cpu=60000

source ~/Data/conda_install/bin/activate ldsc

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/BC/cp_bca_ldsc.txt \
# 	--N 1135754 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/BC/cp_bca_merge

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/BC/bca_ldsc.txt \
# 	--N 1136375 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/BC/bca_merge

# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_bca_merge.sumstats.gz,bca_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.54 \
# 	--pop-prev nan,0.13 \
	--out ~/Data/manansa2/lv_encore/BC/bc_comp_pheno_ldsc

##############################################################################

## ASTHMA

## https://genepi.qimr.edu.au/staff/manuelf/gwas_results/main.html
### 26,582 cases and 300,671controls

#!/usr/bin/env bash

### 

#SBATCH --job-name=ldsc_munge
#SBATCH --mem-per-cpu=60000

source ~/Data/conda_install/bin/activate ldsc

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/ASTHMA/cp_ast_ldsc.txt \
# 	--N 1171290 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/ASTHMA/cp_ast_merge

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/ASTHMA/ast_ldsc.txt \
# 	--N 1171276 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/ASTHMA/ast_merge
# 
# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_ast_merge.sumstats.gz,ast_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--samp-prev nan,0.08 \
# 	--pop-prev nan,0.08 \
	--out ~/Data/manansa2/lv_encore/ASTHMA/ast_comp_pheno_ldsc

##############################################################################

## ALZHEIMERS

## https://www.nature.com/articles/s41588-019-0358-2

#!/usr/bin/env bash

### 

#SBATCH --job-name=alz_ldsc_munge
#SBATCH --mem-per-cpu=60000

source ~/Data/conda_install/bin/activate ldsc

# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/ALZ/cp_alz_ldsc.txt \
# 	--N 1200998 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/ALZ/cp_alz_merge
# 
# ~/Data/manansa2/ldscore/ldsc/munge_sumstats.py \
# 	--sumstats ~/Data/manansa2/lv_encore/ALZ/alz_ldsc.txt \
# 	--N 1200902 \
# 	--merge-alleles ~/Data/manansa2/ldscore/w_hm3.snplist \
# 	--out ~/Data/manansa2/lv_encore/ALZ/alz_merge

# ~/Data/manansa2/ldscore/ldsc/ldsc.py \
# 	--rg cp_alz_merge.sumstats.gz,alz_merge.sumstats.gz \
# 	--ref-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--w-ld-chr  ~/Data/manansa2/ldscore/eur_w_ld_chr/ \
# 	--out ~/Data/manansa2/lv_encore/ALZ/alz_comp_pheno_ldsc
# 	--samp-prev nan,### \
# 	--pop-prev nan,### \
