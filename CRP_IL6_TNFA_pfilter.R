##############################################################################
## Title: P-Value Filtering 
## Purpose: Create phenotype file for use in Encore
## Author: Regina Manansala
## Date Created: 22-August-2019
## Date Modified: 5-September-2019
##############################################################################

library(data.table)
library(dplyr)

## Subset encore index SNPs by pvalue
crp_encore_500kb <- fread("crp_encore_500kb_excl_FIXED.txt")
pval_sub <- crp_encore_500kb[crp_encore_500kb$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())

il6_encore_500kb <- fread("il6_encore_500kb_excl_FIXED.txt")
pval_sub <- il6_encore_500kb[il6_encore_500kb$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())

tnfa_encore_500kb <- fread("tnfa_encore_500kb_excl.txt")
pval_sub <- tnfa_encore_500kb[tnfa_encore_500kb$p.value < 5*10^-8, ] %>% mutate(. ,id = row_number())

## Get index SNPs based on pvalue
pos_sub <- pval_sub
for(i in 1:nrow(pos_sub)){
   pos <- pos_sub[i, "POS"]
   chr <- pos_sub[i, "CHR"]
   min_pos <- pos - 250000
   max_pos <- pos + 250000
   foo <- subset(pos_sub, CHR == chr & POS %between% c(min_pos,max_pos))
   drop_var <- subset(foo, p.value != min(p.value) | duplicated(p.value))
   pos_sub <- pos_sub[!(pos_sub$id %in% drop_var$id),]
   if(i == nrow(pos_sub)) break
}

## Export
write.table(pos_sub[,1:14], "crp_pval_filt_2.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(pos_sub[,1:14], "il6_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)
write.table(pos_sub[,1:14], "tnfa_pval_filt.txt", sep = "\t", row.names = F, col.names = F, quote = F)




// 
// drop_var <- data.frame()
// 
// for(i in 1:nrow(pval_sub)){
//    pos <- pval_sub[i, "POS"]
//    chr <- pval_sub[i, "CHR"]
//    min_pos <- pos - 250000
//    max_pos <- pos + 250000
//    pos_sub <- subset(pval_sub, CHR == chr & POS %between% c(min_pos,max_pos))
//    drop_var <- subset(pos_sub, p.value != min(p.value) | duplicated(p.value)) %>% rbind(., drop_var)
//    #pos_sub <- pos_sub[!(pos_sub$id %in% drop_var$id),]
//    if(i == nrow(pos_sub)) break
// }
