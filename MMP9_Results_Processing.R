##############################################################################
## Title: MMP9 Encore Result Filtering
## Purpose: Subset index SNPS + SNPs within 500kb (+/- 250kb) and filter 
##		top 50 hits
## Author: Regina Manansala
## Date Created: 24-October-2019
## Date Modified: 24-October-2019
##############################################################################

library(data.table)
library(dplyr)
library(qqman)

## Read encore results
mmp9_encore <- fread("mmp9_results.txt")

## Read reference GWAS results
mmp9_ref <- fread("gwas-association-downloaded_2019-09-20-EFO_0004746-LPPLA2.tsv") %>%
	subset(., (`P-VALUE (TEXT)` %in% c("(mass)", "(Mass concentrations)"))) %>%
	select(., c("SNPS", "CHR_ID", "CHR_POS"))
pub_snp <- rbind(c("rs34159425", 6, 46705176), c("rs1051931", 6, 46705206), c("rs76863441", 6, 46709361), c("rs144983904", 6, 46716515), c("rs140020965", 6, 46709337), c("rs142974898", 6, 46722781)) %>% 
	data.frame() %>%
	mutate(., X2 = as.character(X2) %>% as.integer(X2), X3 = as.character(X3) %>% as.integer(X3))
names(pub_snp) <- names(mmp9_ref)
mmp9_ref <- rbind(mmp9_ref, pub_snp) #%>% mutate(., CHR_ID = as.numeric(CHR_ID), CHR_POS = as.numeric(CHR_POS))

## Create new Chromosome column without "chr" prefix
mmp9_encore$CHR_2 <- sub("chr", "", mmp9_encore$CHR)
mmp9_encore$CHR_2[mmp9_encore$CHR_2 == "X"] <- "23"
mmp9_encore$CHR_2 <- as.numeric(mmp9_encore$CHR_2)

## Join encore and reference results to get index SNPs
mmp9_sub <- left_join(mmp9_ref[, c("CHR_ID", "CHR_POS", "SNPS")], mmp9_encore, by=c("CHR_ID"="CHR_2", "CHR_POS"="POS"))
mmp9_sub <- subset(mmp9_sub, !is.na(p.value))
write.table(mmp9_sub, "mmp9_gwascat.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Get list of chr and pos from index SNPs
chr <- as.list(mmp9_sub$CHR_ID) %>% unlist()
pos <- as.list(mmp9_sub$CHR_POS) %>% unlist()

# mmp9_pos_excl <- mmp9_encore
# n_excl <- 0
# for(i in 1:44){
# 	min_pos <- pos[i] - 250000
# 	max_pos <- pos[i] + 250000
# 	start_n <- nrow(mmp9_pos_excl)
# 	mmp9_pos_excl <- subset(mmp9_pos_excl, !(CHR_2 == chr[i] & POS > min_pos & POS < max_pos))
# 	end_n <- nrow(mmp9_pos_excl)
# 	change_n <- start_n - end_n
# 	n_excl <- n_excl + change_n
# }

## Subset encore results to index SNPs + SNPs within 500kb (+/- 250kb)
mmp9_pos_excl <- mmp9_encore
n_excl <- 0
for(i in 1:5){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	n_excl <- n_excl + nrow(mmp9_pos_excl[mmp9_pos_excl$POS > min_pos & mmp9_pos_excl$POS < max_pos & mmp9_pos_excl$CHR_2 == chr[i],])
	print(nrow(mmp9_pos_excl[mmp9_pos_excl$POS > min_pos & mmp9_pos_excl$POS < max_pos & mmp9_pos_excl$CHR_2 == chr[i],]))
	mmp9_pos_excl <- subset(mmp9_pos_excl, !(CHR_2 == chr[i] & POS > min_pos & POS < max_pos))
}

# n_excl <- 0
# for(i in 1:44){
# 	min_pos <- pos[i] - 250000
# 	max_pos <- pos[i] + 250000
# 	n_excl <- n_excl + nrow(mmp9_encore[mmp9_encore$POS > min_pos & mmp9_encore$POS < max_pos & mmp9_encore$CHR_2 == chr[i],])
# }

## Export results
write.table(mmp9_pos_excl, "mmp9_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Export plots of pvalue distributions
jpeg("mmp9_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(mmp9_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("mmp9_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(mmp9_pos_excl$p.value)
dev.off()

## Filter p-values
mmp9_pval_sub <- mmp9_pos_excl[mmp9_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
for(i in 1:nrow(mmp9_pval_sub)){
   pos <- mmp9_pval_sub[i, "POS"]
   chr <- mmp9_pval_sub[i, "CHR"]
   min_pos <- pos - 250000
   max_pos <- pos + 250000
   foo <- subset(mmp9_pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, p.value != min(p.value) | duplicated(p.value))
   mmp9_pval_sub <- mmp9_pval_sub[!(mmp9_pval_sub$id %in% drop_var$id),]
   if(i == nrow(mmp9_pval_sub)) break
}
write.table(mmp9_pval_sub[,1:14], "mmp9_mass_pval_filt.txt", sep = "\t", row.names = F, col.names = T, quote = F)
