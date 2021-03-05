##############################################################################
## Title: Lp-PLA2 Activity & Mass Encore Result Filtering
## Purpose: Subset index SNPS + SNPs within 500kb (+/- 250kb) and filter 
##		top 50 hits
## Author: Regina Manansala
## Date Created: 26-September-2019
## Date Modified: 22-March-2020
##############################################################################

library(data.table)
library(dplyr)
library(qqman)

## Read encore results
lppla2_encore <- fread("lppla2_act_wFHS_results.txt")
lppla2_encore <- fread("lppla2_mass_wFHS_results.txt")

## Create new Chromosome column without "chr" prefix
lppla2_encore$CHR_2 <- sub("chr", "", lppla2_encore$CHR)
lppla2_encore$CHR_2[lppla2_encore$CHR_2 == "X"] <- "23"
lppla2_encore$CHR_2 <- as.numeric(lppla2_encore$CHR_2)

## Read reference GWAS results
lppla2_ref <- fread("../lppla2_act/gwas-association-downloaded_2019-09-20-EFO_0004746-LPPLA2.tsv") %>%
	subset(., !(`P-VALUE (TEXT)` %in% c("(mass)", "(Mass concentrations)"))) %>%
	select(., c("SNPS", "CHR_ID", "CHR_POS"))
pub_snp <- rbind(c("rs34159425", 6, 46705176), c("rs1051931", 6, 46705206), c("rs76863441", 6, 46709361), c("rs144983904", 6, 46716515), c("rs140020965", 6, 46709337), c("rs142974898", 6, 46722781)) %>% 
	data.frame() %>%
	mutate(., X2 = as.character(X2) %>% as.integer(X2), X3 = as.character(X3) %>% as.integer(X3))

lppla2_ref <- fread("../lppla2_act/gwas-association-downloaded_2019-09-20-EFO_0004746-LPPLA2.tsv") %>%
	subset(., (`P-VALUE (TEXT)` %in% c("(mass)", "(Mass concentrations)"))) %>%
	select(., c("SNPS", "CHR_ID", "CHR_POS"))

names(pub_snp) <- names(lppla2_ref)
lppla2_ref <- rbind(lppla2_ref, pub_snp) #%>% mutate(., CHR_ID = as.numeric(CHR_ID), CHR_POS = as.numeric(CHR_POS))

## Join encore and reference results to get index SNPs
lppla2_sub <- left_join(lppla2_ref[, c("CHR_ID", "CHR_POS", "SNPS")], lppla2_encore, by=c("CHR_ID"="CHR_2", "CHR_POS"="POS"))
lppla2_sub <- subset(lppla2_sub, !is.na(p.value))
write.table(lppla2_sub, "lppla2_act_topmed_032020.txt", sep = "\t", row.names = F, col.names=T, quote = F)
write.table(lppla2_sub, "lppla2_mass_topmed_032020.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Get list of chr and pos from index SNPs
chr <- as.list(lppla2_sub$CHR_ID) %>% unlist()
pos <- as.list(lppla2_sub$CHR_POS) %>% unlist()

## Subset encore results to index SNPs + SNPs within 500kb (+/- 250kb)
# lppla2_pos_excl <- lppla2_encore
# n_excl <- 0
# for(i in 1:44){
# 	min_pos <- pos[i] - 250000
# 	max_pos <- pos[i] + 250000
# 	start_n <- nrow(lppla2_pos_excl)
# 	lppla2_pos_excl <- subset(lppla2_pos_excl, !(CHR_2 == chr[i] & POS > min_pos & POS < max_pos))
# 	end_n <- nrow(lppla2_pos_excl)
# 	change_n <- start_n - end_n
# 	n_excl <- n_excl + change_n
# }

lppla2_pos_excl <- lppla2_encore
n_excl <- 0
for(i in 1:41){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	n_excl <- n_excl + nrow(lppla2_pos_excl[lppla2_pos_excl$POS > min_pos & lppla2_pos_excl$POS < max_pos & lppla2_pos_excl$CHR_2 == chr[i],])
	print(nrow(lppla2_pos_excl[lppla2_pos_excl$POS > min_pos & lppla2_pos_excl$POS < max_pos & lppla2_pos_excl$CHR_2 == chr[i],]))
	lppla2_pos_excl <- subset(lppla2_pos_excl, !(CHR_2 == chr[i] & POS > min_pos & POS < max_pos))
}

# n_excl <- 0
# for(i in 1:44){
# 	min_pos <- pos[i] - 250000
# 	max_pos <- pos[i] + 250000
# 	n_excl <- n_excl + nrow(lppla2_encore[lppla2_encore$POS > min_pos & lppla2_encore$POS < max_pos & lppla2_encore$CHR_2 == chr[i],])
# }

## Export results
write.table(lppla2_pos_excl, "lppla2_act_500kb_excl_032020.txt", sep = "\t", row.names = F, col.names=T, quote = F)
write.table(lppla2_pos_excl, "lppla2_mass_500kb_excl_032020.txt", sep = "\t", row.names = F, col.names=T, quote = F)

## Export plots of pvalue distributions
jpeg("lppla2_act_excl_manhattan_032020.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(lppla2_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("lppla2_act_excl_qq_032020.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(lppla2_pos_excl$p.value)
dev.off()

jpeg("lppla2_mass_excl_manhattan_032020.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(lppla2_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("lppla2_mass_excl_qq_032020.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(lppla2_pos_excl$p.value)
dev.off()

## Filter p-values
lppla2_pval_sub <- lppla2_pos_excl[lppla2_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
for(i in 1:nrow(lppla2_pval_sub)){
   pos <- lppla2_pval_sub[i, "POS"]
   chr <- lppla2_pval_sub[i, "CHR"]
   min_pos <- pos - 250000
   max_pos <- pos + 250000
   foo <- subset(lppla2_pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, p.value != min(p.value) | duplicated(p.value))
   lppla2_pval_sub <- lppla2_pval_sub[!(lppla2_pval_sub$id %in% drop_var$id),]
   if(i == nrow(lppla2_pval_sub)) break
}
write.table(lppla2_pval_sub[,1:14], "lppla2_act_pval_filt_032020.txt", sep = "\t", row.names = F, col.names = T, quote = F)



