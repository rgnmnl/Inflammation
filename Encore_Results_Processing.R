##############################################################################
## Title: Encore Result Filtering
## Purpose: Subset index SNPS + SNPs within 500kb (+/- 250kb) and filter 
##		top 50 hits
## Author: Regina Manansala
## Date Created: 15-July-2019
## Date Modified: 30-September-2020
##############################################################################

library(data.table)
library(dplyr)
library(qqman)

## Read CRP Encore results and reference GWAS results
crp_encore_full <- fread("results.txt")
crp_ref <- fread("crp_ref_results.txt")

## Create new Chromosome column without "chr" prefix
crp_encore_full$CHR_2 <- sub("chr", "", crp_encore_full$CHR)
crp_encore_full$CHR_2[crp_encore_full$CHR_2 == "X"] <- "23"
crp_encore_full$CHR_2 <- as.numeric(crp_encore_full$CHR_2)

## Join encore and reference results to get index SNPs
crp_sub <- left_join(crp_ref_results[, c("Chr", "Position_hg38", "Variant", "Coded_Allele")], crp_encore_full, by=c("Chr"="CHR_2", "Position_hg38"="POS"))
crp_sub_2 <- subset(crp_sub, !is.na(p.value))
crp_sub_2[, c("Chr", "Position_hg38", "Variant", "Coded_Allele", "CHR", "Allele1", "Allele2")]

## Get list of chr and pos from index SNPs
chr <- as.list(crp_sub_2$Chr) %>% unlist()
pos <- as.list(crp_sub_2$Position_hg38) %>% unlist()

## Subset encore results to index SNPs + SNPs within 500kb (+/- 250kb)
pos_excl <- crp_encore_full
for(i in 1:57){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	pos_excl <- subset(pos_excl, !(CHR == chr[i] &POS > min_pos & POS < max_pos))
}

## Check number of rows removed from encore results
n_excl <- 0
for(i in 1:57){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	n_excl <- n_excl + nrow(crp_encore[crp_encore$POS > min_pos & crp_encore$POS < max_pos & crp_encore$CHR == chr[i],])
}

## Checks
dim(pos_excl[pos_excl$POS > min_pos & pos_excl$POS < max_pos & pos_excl$CHR_2 == chr[i],])
dim(crp_encore_full[crp_encore_full$POS > min_pos & crp_encore_full$POS < max_pos & crp_encore_full$CHR_2 == chr[i],])

## Export plots of pvalue distributions
jpeg("crp_500kb_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("crp_500kb_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(pos_excl$p.value)
dev.off()

## Export results
write.table(pos_excl, "crp_encore_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Get top 50 hits from Encore and export
top_50_encore <- pos_excl[order(pos_excl$p.value),]
top_50_encore <- top_50_encore[1:50,]
write.table(top_50_encore, "crp_encore_top50.txt", sep = "\t", row.names = F, col.names=T, quote = F)

##############################################################################
############################### TNFA Filtering ###############################
##############################################################################

## Read CRP Encore results and reference GWAS results
tnfa_encore_full <- fread("results.txt")
tnfa_ref <- fread("gwas-association-TNFA.tsv")

## Get list of chr and pos from index SNPs
chr <- as.list(tnfa_ref$CHR_ID) %>% unlist()
pos <- as.list(tnfa_ref$CHR_POS) %>% unlist()

## Subset encore results to index SNPs + SNPs within 500kb (+/- 250kb)
pos_excl <- tnfa_encore_full
for(i in 1:32){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	pos_excl <- subset(pos_excl, !(CHR == chr[i] &POS > min_pos & POS < max_pos))
}

## Check number of rows removed from encore results
n_excl <- 0
for(i in 1:32){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	n_excl <- n_excl + nrow(tnfa_encore[tnfa_encore$POS > min_pos & tnfa_encore$POS < max_pos & tnfa_encore$CHR == chr[i],])
}

## Export plots of pvalue distributions
jpeg("tnfa_500kb_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("tnfa_500kb_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(pos_excl$p.value)
dev.off()

##############################################################################
########################  CRP RACE STRATIFIED RESULTS ########################
##############################################################################

## Read CRP race stratified results from encore
crp_white <- fread("crp_white_results.txt")
crp_black <- fread("crp_black_results.txt")
crp_hisp <- fread("crp_hisp_results.txt")

## Create new Chromosome column without "chr" prefix
crp_white$CHR_2 <- sub("chr", "", crp_white$CHR)
crp_white$CHR_2[crp_white$CHR_2 == "X"] <- "23"
crp_white$CHR_2 <- as.numeric(crp_white$CHR_2)

crp_black$CHR_2 <- sub("chr", "", crp_black$CHR)
crp_black$CHR_2[crp_black$CHR_2 == "X"] <- "23"
crp_black$CHR_2 <- as.numeric(crp_black$CHR_2)

crp_hisp$CHR_2 <- sub("chr", "", crp_hisp$CHR)
crp_hisp$CHR_2[crp_hisp$CHR_2 == "X"] <- "23"
crp_hisp$CHR_2 <- as.numeric(crp_hisp$CHR_2)

## Join encore and reference results to get index SNPs
crp_white_sub <- left_join(crp_ref[, c("Chr", "Position_hg38", "Variant", "Coded_Allele")], crp_white, by=c("Chr"="CHR_2", "Position_hg38"="POS"))
crp_white_sub_2 <- subset(crp_white_sub, !is.na(p.value))
write.table(crp_white_sub_2, "crp_white_topmed.txt", sep = "\t", row.names = F, col.names=T, quote = F)

crp_black_sub <- left_join(crp_ref[, c("Chr", "Position_hg38", "Variant", "Coded_Allele")], crp_black, by=c("Chr"="CHR_2", "Position_hg38"="POS"))
crp_black_sub_2 <- subset(crp_black_sub, !is.na(p.value))
write.table(crp_black_sub_2, "crp_black_topmed.txt", sep = "\t", row.names = F, col.names=T, quote = F)

crp_hisp_sub <- left_join(crp_ref[, c("Chr", "Position_hg38", "Variant", "Coded_Allele")], crp_hisp, by=c("Chr"="CHR_2", "Position_hg38"="POS"))
crp_hisp_sub_2 <- subset(crp_hisp_sub, !is.na(p.value))
write.table(crp_hisp_sub_2, "crp_hisp_topmed.txt", sep = "\t", row.names = F, col.names=T, quote = F)

# crp_white_sub_2[, c("Chr", "Position_hg38", "Variant", "Coded_Allele", "CHR", "Allele1", "Allele2")]

## Get list of chr and pos from index SNPs
chr <- as.list(crp_white_sub_2$Chr) %>% unlist()
pos <- as.list(crp_white_sub_2$Position_hg38) %>% unlist()

## Subset encore results to index SNPs + SNPs within 500kb (+/- 250kb)
white_pos_excl <- crp_white
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	white_pos_excl <- subset(white_pos_excl, !(CHR_2 == chr[i] &POS > min_pos & POS < max_pos))
}

## Check number of rows removed from encore results
n_excl <- 0
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	n_excl <- n_excl + nrow(crp_white[crp_white$POS > min_pos & crp_white$POS < max_pos & crp_white$CHR_2 == chr[i],])
}

write.table(white_pos_excl, "crp_white_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Export plots of pvalue distributions
jpeg("crp_white_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(white_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("crp_white_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(white_pos_excl$p.value)
dev.off()

## Subset encore results to index SNPs + SNPs within 500kb (+/- 250kb)
black_pos_excl <- crp_black
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	black_pos_excl <- subset(black_pos_excl, !(CHR_2 == chr[i] &POS > min_pos & POS < max_pos))
}
write.table(black_pos_excl, "crp_black_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Export plots of pvalue distributions
jpeg("crp_black_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(black_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("crp_black_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(black_pos_excl$p.value)
dev.off()

## Subset encore results to index SNPs + SNPs within 500kb (+/- 250kb)
hisp_pos_excl <- crp_hisp
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	hisp_pos_excl <- subset(hisp_pos_excl, !(CHR_2 == chr[i] &POS > min_pos & POS < max_pos))
}
write.table(hisp_pos_excl, "crp_hisp_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Export plots of pvalue distributions
jpeg("crp_hisp_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(hisp_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("crp_hisp_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(hisp_pos_excl$p.value)
dev.off()

## Filter p-values
white_pval_sub <- white_pos_excl[white_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
for(i in 1:nrow(white_pval_sub)){
   pos <- white_pval_sub[i, "POS"]
   chr <- white_pval_sub[i, "CHR"]
   min_pos <- pos - 250000
   max_pos <- pos + 250000
   foo <- subset(white_pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, p.value != min(p.value) | duplicated(p.value))
   white_pval_sub <- white_pval_sub[!(white_pval_sub$id %in% drop_var$id),]
   if(i == nrow(white_pval_sub)) break
}
write.table(white_pval_sub[,1:14], "crp_white_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)

## Filter p-valuesblack_pval_sub <- black_pos_excl[black_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
for(i in 1:nrow(black_pval_sub)){
   pos <- black_pval_sub[i, "POS"]
   chr <- black_pval_sub[i, "CHR"]
   min_pos <- pos - 250000
   max_pos <- pos + 250000
   foo <- subset(black_pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, p.value != min(p.value) | duplicated(p.value))
   black_pval_sub <- black_pval_sub[!(black_pval_sub$id %in% drop_var$id),]
   if(i == nrow(black_pval_sub)) break
}
write.table(black_pval_sub[,1:14], "crp_black_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)


## Filter p-valueshisp_pval_sub <- hisp_pos_excl[hisp_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
for(i in 1:nrow(hisp_pval_sub)){
   pos <- hisp_pval_sub[i, "POS"]
   chr <- hisp_pval_sub[i, "CHR"]
   min_pos <- pos - 250000
   max_pos <- pos + 250000
   foo <- subset(hisp_pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, p.value != min(p.value) | duplicated(p.value))
   hisp_pval_sub <- hisp_pval_sub[!(hisp_pval_sub$id %in% drop_var$id),]
   if(i == nrow(hisp_pval_sub)) break
}
write.table(hisp_pval_sub[,1:14], "crp_hisp_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)


##############

