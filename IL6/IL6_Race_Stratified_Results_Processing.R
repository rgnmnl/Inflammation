##############################################################################
## Title: IL6 Encore Result Filtering
## Purpose: Subset index SNPS + SNPs within 500kb (+/- 250kb) and filter 
##		top 50 hits for IL6 stratified gwas results
## Author: Regina Manansala
## Date Created: 30-September-2019
## Date Modified: 08-October-2019
##############################################################################

library(data.table)
library(dplyr)
library(qqman)

il6_gwas <- fread("il6_gwas_snps_complete.txt")

il6_whi <- fread("Stratified_Analysis/IL6_WHI_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_mesa <- fread("Stratified_Analysis/IL6_MESA_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_fhs <- fread("Stratified_Analysis/IL6_FHS_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_fhs1 <- fread("Stratified_Analysis/IL6_FHS1_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_fhs3 <- fread("Stratified_Analysis/IL6_FHS3_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_copd <- fread("Stratified_Analysis/IL6_COPDGene_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_chs <- fread("Stratified_Analysis/IL6_CHS_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_cfs <- fread("Stratified_Analysis/IL6_CFS_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_cfs <- fread("Stratified_Analysis/IL6_CFS_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_cardia <- fread("Stratified_Analysis/IL6_CARDIA_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453) 
il6_whi_elisa <- fread("Stratified_Analysis/IL6_WHI_ELISA_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)
il6_whi_hs <- fread("Stratified_Analysis/IL6_WHI_HS_results.txt") %>% subset(., CHR == "chr2" & POS == 113083453)

il6_mesa$STUDY <- "MESA"
il6_fhs$STUDY <- "FHS"
il6_fhs1$STUDY <- "FHS - Cohort 1"
il6_fhs3$STUDY <- "FHS - Cohort 3"
il6_copd$STUDY <- "COPDGene"
il6_chs$STUDY <- "CHS"
il6_cfs$STUDY <- "CFS"
il6_cardia$STUDY <- "CARDIA"
il6_whi_elisa$STUDY <- "WHI - HS ELISA"
il6_whi_hs$STUDY <- "WHI - HS ASSAY"
il6_whi$STUDY <- "WHI"

il6_rs6734238 <- rbind(il6_mesa, il6_fhs, il6_fhs1, il6_fhs3, il6_copd, il6_chs, il6_cfs, il6_cardia[, -8], il6_whi_elisa, il6_whi_hs)

write.table(il6_rs6734238, "Stratified_Analysis/il6_rs6734238.txt", sep = "\t", row.names = F, col.names = T, quote = F)

###################

il6_white <- fread("Stratified_Analysis/IL6_white_results.txt") #%>% subset(., CHR == "chr2" & POS == 113083453)
il6_black <- fread("Stratified_Analysis/IL6_BLACK_results.txt") #%>% subset(., CHR == "chr2" & POS == 113083453)
il6_hisp <- fread("Stratified_Analysis/IL6_HISPANIC_results.txt") #%>% subset(., CHR == "chr2" & POS == 113083453)

il6_ref <- fread("il6_gwas_snps_complete.txt")

il6_white_chr2 <- subset(il6_white, CHR == 'chr2' & POS %between% c(min,max))
il6_black_chr2 <- subset(il6_black, CHR == 'chr2' & POS %between% c(min,max))
il6_hisp_chr2 <- subset(il6_hisp, CHR == 'chr2' & POS %between% c(min,max))

write.table(il6_white_chr2, "Stratified_Analysis/il6_white_500kb_rsID.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(il6_black_chr2, "Stratified_Analysis/il6_black_500kb_rsID.txt", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(il6_hisp_chr2, "Stratified_Analysis/il6_hisp_500kb_rsID.txt", sep = "\t", row.names = F, col.names = T, quote = F)


jpeg("il6_white_chr2_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
#manhattan(white_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
manhattan(subset(il6_white, CHR_2 == 2), chr="CHR_2", bp="POS", p="p.value", cex = 0.1, main = "CHR 2")
dev.off()

// il6_white$race <- "white"
// il6_black$race <- "BLACK"
// il6_hispanic$race <- "HISPANIC"
// 
// il6_rs6734238_race <- rbind(il6_white, il6_black, il6_hispanic)
// 
// write.table(il6_rs6734238_race, "Stratified_Analysis/il6_rs6734238_race.txt", sep = "\t", row.names = F, col.names = T, quote = F)

il6_white$CHR_2 <- sub("chr", "", il6_white$CHR)
il6_white$CHR_2[il6_white$CHR_2 == "X"] <- "23"
il6_white$CHR_2 <- as.numeric(il6_white$CHR_2)

il6_black$CHR_2 <- sub("chr", "", il6_black$CHR)
il6_black$CHR_2[il6_black$CHR_2 == "X"] <- "23"
il6_black$CHR_2 <- as.numeric(il6_black$CHR_2)

il6_hisp$CHR_2 <- sub("chr", "", il6_hisp$CHR)
il6_hisp$CHR_2[il6_hisp$CHR_2 == "X"] <- "23"
il6_hisp$CHR_2 <- as.numeric(il6_hisp$CHR_2)

il6_white_sub <- left_join(il6_ref, il6_white, by=c("chr_2"="CHR", "pos"="POS"))
il6_white_sub_2 <- subset(il6_white_sub, !is.na(p.value))
write.table(il6_white_sub_2, "il6_white_gwascat.txt", sep = "\t", row.names = F, col.names=T, quote = F)

il6_black_sub <- left_join(il6_ref, il6_black, by=c("chr_2"="CHR", "pos"="POS"))
il6_black_sub_2 <- subset(il6_black_sub, !is.na(p.value))
write.table(il6_black_sub_2, "il6_black_gwascat.txt", sep = "\t", row.names = F, col.names=T, quote = F)

il6_hisp_sub <- left_join(il6_ref, il6_hisp, by=c("chr_2"="CHR", "pos"="POS"))
il6_hisp_sub_2 <- subset(il6_hisp_sub, !is.na(p.value))
write.table(il6_hisp_sub_2, "il6_hisp_gwascat.txt", sep = "\t", row.names = F, col.names=T, quote = F)

##

chr <- as.list(il6_ref$chr) %>% unlist()
pos <- as.list(il6_pos$pos) %>% unlist()

white_pos_excl <- il6_white
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	white_pos_excl <- subset(white_pos_excl, !(CHR_2 == chr[i] &POS > min_pos & POS < max_pos))
}

n_excl <- 0
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	n_excl <- n_excl + nrow(il6_white[il6_white$POS > min_pos & il6_white$POS < max_pos & il6_white$CHR_2 == chr[i],])
}

write.table(white_pos_excl, "il6_white_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

jpeg("il6_white_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(white_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("il6_white_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(white_pos_excl$p.value)
dev.off()

black_pos_excl <- il6_black
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	black_pos_excl <- subset(black_pos_excl, !(CHR_2 == chr[i] &POS > min_pos & POS < max_pos))
}
write.table(black_pos_excl, "il6_black_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

jpeg("il6_black_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(black_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("il6_black_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(black_pos_excl$p.value)
dev.off()

hisp_pos_excl <- il6_hisp
for(i in 1:55){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	hisp_pos_excl <- subset(hisp_pos_excl, !(CHR_2 == chr[i] &POS > min_pos & POS < max_pos))
}
write.table(hisp_pos_excl, "il6_hisp_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

jpeg("il6_hisp_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(hisp_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("il6_hisp_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(hisp_pos_excl$p.value)
dev.off()

##

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
write.table(white_pval_sub[,1:14], "il6_white_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)

black_pval_sub <- black_pos_excl[black_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
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
write.table(black_pval_sub[,1:14], "il6_black_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)


hisp_pval_sub <- hisp_pos_excl[hisp_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
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
write.table(hisp_pval_sub[,1:14], "il6_hisp_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)


################

il6_cond <- fread("Stratified_Analysis/IL6_CONDITIONAL_results.txt")

il6_cond$CHR_2 <- sub("chr", "", il6_cond$CHR)
il6_cond$CHR_2[il6_cond$CHR_2 == "X"] <- "23"
il6_cond$CHR_2 <- as.numeric(il6_cond$CHR_2)

il6_cond_sub <- left_join(il6_ref, il6_cond, by=c("chr_2"="CHR", "pos"="POS")) %>% subset(!is.na(p.value))
write.table(il6_cond_sub, "il6_cond_gwascat.txt", sep = "\t", row.names = F, col.names=T, quote = F)

chr <- as.list(il6_ref$chr) %>% unlist()
pos <- as.list(il6_ref$pos) %>% unlist()

cond_pos_excl <- il6_cond
for(i in 1:84){
	min_pos <- pos[i] - 250000
	max_pos <- pos[i] + 250000
	cond_pos_excl <- subset(cond_pos_excl, !(CHR_2 == chr[i] &POS > min_pos & POS < max_pos))
}

write.table(cond_pos_excl, "il6_cond_500kb_excl.txt", sep = "\t", row.names = F, col.names=T, quote = F)

jpeg("il6_cond_excl_manhattan.jpg", height=7, width=10, units='in', quality=100, res=500)
manhattan(cond_pos_excl, chr="CHR_2", bp="POS", p="p.value", cex = 0.6)
dev.off()

jpeg("il6_cond_excl_qq.jpg", height=10, width=10, units='in', quality=100, res=500)
qq(cond_pos_excl$p.value)
dev.off()

cond_pval_sub <- cond_pos_excl[cond_pos_excl$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())
for(i in 1:nrow(cond_pval_sub)){
   pos <- cond_pval_sub[i, "POS"]
   chr <- cond_pval_sub[i, "CHR"]
   min_pos <- pos - 250000
   max_pos <- pos + 250000
   foo <- subset(cond_pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, p.value != min(p.value) | duplicated(p.value))
   cond_pval_sub <- cond_pval_sub[!(cond_pval_sub$id %in% drop_var$id),]
   if(i == nrow(cond_pval_sub)) break
}

write.table(cond_pval_sub[,1:14], "il6_cond_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)

min <- 113083453 - 250000
max <- 113083453 + 250000
il6_cond_chr2 <- subset(il6_cond, CHR == 'chr2' & POS %between% c(min,max))
write.table(il6_cond_chr2, "Stratified_Analysis/il6_cond_500kb_rsID.txt", sep = "\t", row.names = F, col.names = T, quote = F)
