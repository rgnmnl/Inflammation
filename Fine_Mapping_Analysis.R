##############################################################################
## Title: Fine Mapping Analysis
## Version: 2
## Reference: https://www.biorxiv.org/content/10.1101/2020.01.17.910497v1
## Author: Regina Manansala
## Date Created: 09-October-2019
## Date Modified: 14-October-2019
##############################################################################

library(data.table)
library(dplyr)

bayesfactor <- function(file_path, chr, pos) {
  min <- pos - 250000
  max <- pos + 250000
  encore_result <- fread(file_path) %>% subset(., CHR == chr & POS %between% c(min,max))
  
  ## Calculate variance
  encore_result$V <- (encore_result$SE)^2
  ## Calculate Z score
  encore_result$Z <- encore_result$BETA/encore_result$SE
  ## Specify prior
  W <- .04
  ## Calculate approximate Bayes Factor
  encore_result$aBF <- sqrt(encore_result$V/encore_result$V+W) * exp((encore_result$Z^2 * W)/(2*(encore_result$V + W)))
  ## Calculate PPi. Sort and subset data
  encore_result$PPi <- encore_result$aBF/sum(encore_result$aBF)
  
  encore_result <- encore_result[order(-encore_result$PPi),]
  
  N <- 0
  keep_rows <- data.frame()
  for(i in 1:nrow(encore_result)){
    PPi <- encore_result[i]$PPi
    N <- N + PPi
    keep_rows <- rbind(keep_rows, encore_result[i])
    if(N > .95) break
  }
  return(keep_rows)
}

