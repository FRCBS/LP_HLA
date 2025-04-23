
library(tidyverse)
library(data.table)
library(glue)
library(corrplot)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(knitr)
library(scales)
library(egg)
library(readxl)


## Functions

# make a new covar table with selected HLA alleles
# expects covarites in 'covar' and hla alleles in 'hla' variables
# adds hla alleles as new covariates
makeNewCov <- function(x, alleles) {
  x <- inner_join(covars, hla[, c('IID', alleles)])
  colnames(x) <- gsub('*', '.', colnames(x), fixed=T) %>% gsub(':', '', ., fixed=T)
  return(x)
}


# create plink formatted genotype file
# input matrix d has IID sample ID col, and the rest are genotype dosages with descriptive col names
createPlinkGeno <- function(d, outname) {
  
  iid <- d$IID
  d <- d %>% dplyr::select(-IID)
  
  d <- d %>% mutate_all(as.character)
  
  d[d=='0'] <- 'absent absent'
  d[d=='1'] <- 'absent present'
  d[d=='2'] <- 'present present'
  
  haplos.ped <- data.frame(iid,
                           iid,
                           0,0,0,0,
                           d
  )
  haplos.map <- data.frame(rep('6', ncol(d)),
                           colnames(d),
                           rep(0, ncol(d)),
                           rep(32123123, ncol(d)) # the position is just something in the MHC
  )
  
  
  fwrite(haplos.ped, glue('data/genotypes/{outname}_haplos.ped'), sep=' ', col.names=F, quote=F)
  fwrite(haplos.map, glue('data/genotypes/{outname}_haplos.map'), sep=' ', col.names=F, quote=F)
  
  system(glue("plink --file data/genotypes/{outname}_haplos --make-bed --out data/genotypes/{outname}_haplos"))
  #system(glue("plink2 --bfile data/genotypes/{outname}_haplos --freq"))
}


# define homozygotes, heterozygotes and additive genotypes from a HLA allele pair 
# expects HLA allele dosages matrix to be in variable 'hla'
# outputs to 'data/genotypes/output.name'
makeHLAhom <- function(alleles, output.name) {
  
  h1 <- hla[, alleles]
  sum(h1[,1]>1.5 & h1[,2]>1.5)
  sum(h1[,1]>0.5 & h1[,2]>0.5 & h1[,1]<1.5 & h1[,2]<1.5)
  sum(h1[,1]<0.5 & h1[,2]<0.5)
  
  h1.hom <- rep(0, nrow(hla))
  h1.hom[h1[, 1]>1.5 & h1[, 2]>1.5] <- 1 
  
  h1.het <- rep(0, nrow(hla))
  h1.het[h1[,1]>0.5 & h1[,2]>0.5 & h1[,1]<1.5 & h1[,2]<1.5] <- 1 
  
  h1.add <- rep(0, nrow(hla))
  h1.add[h1[,1]>1.5 & h1[,2]>1.5] <- 2
  h1.add[h1[,1]>0.5 & h1[,2]>0.5 & h1[,1]<1.5 & h1[,2]<1.5] <- 1 
  
  h1 <- data.frame(IID=hla$IID, 
                   X_hom = h1.hom,
                   X_het = h1.het,
                   X = h1.add)
  
  colnames(h1)[2:4] <- paste(
    gsub('*', '.', alleles[1], fixed=T) %>% gsub(':', '', ., fixed=T),
    gsub('*', '.', alleles[2], fixed=T) %>% gsub(':', '', ., fixed=T),
    sep='_') %>% 
    paste(., c('hom', 'het', 'add'), sep='_')
  
  createPlinkGeno(h1, output.name)
  
}

# convert log p -values bacto to p
convert.p <- function(x) return(10^-x)

# FDR correction
fdr_by <- function(x) p.adjust(x, method='BY') %>% return


# process regenie assoc results of HLA alleles
process_result <- function(x, p) {
  x <- fread(x) %>% filter(!grepl('DRB3|DRB4|DRB5', ID))
  x <- x %>% mutate(P=convert.p(LOG10P)) %>% 
    mutate(FDR_BY=fdr_by(P)) %>% #filter(FDR_BY<0.01) %>% 
    arrange(1/LOG10P) %>%
    dplyr::select(c(ID, A1FREQ, N, TEST, BETA, SE, LOG10P, P, FDR_BY))
  x$CI95  <- 1.96*x$SE
  x$UPPER <- x$BETA+x$CI95
  x$LOWER <- x$BETA-x$CI95
  x$OR <- exp(x$BETA)
  x$OR_UPPER <- exp(x$UPPER)
  x$OR_LOWER <- exp(x$LOWER)
  x$OR_95CI <- paste0(x$OR_LOWER %>% round(2), ' - ', x$OR_UPPER %>% round(2))
  
  x <- data.frame(x, x$ID %>% str_split_fixed(., '_', 2))
  x$Gene <- str_split_fixed(x$ID, '\\*', 2)[, 1]
  return(data.frame(x, Pheno=p) %>% dplyr::select(-c(X1, X2)))
}


# results processing for DQB1 pair analysis
process_result_dqb1pair <- function(x, p) {
  x <- fread(x)
  x <- x %>% mutate(P=convert.p(LOG10P)) %>% 
    mutate(FDR_BY=fdr_by(P)) %>% #filter(FDR_BY<0.01) %>% 
    arrange(1/LOG10P) %>%
    dplyr::select(c(ID, A1FREQ, N, TEST, BETA, SE, LOG10P, P, FDR_BY))
  x$CI95  <- 1.96*x$SE
  x$UPPER <- x$BETA+x$CI95
  x$LOWER <- x$BETA-x$CI95
  x$OR <- exp(x$BETA)
  x$OR_UPPER <- exp(x$UPPER)
  x$OR_LOWER <- exp(x$LOWER)
  x$OR_95CI <- paste0(x$OR_LOWER %>% round(2), ' - ', x$OR_UPPER %>% round(2))
  
  x <- data.frame(x, x$ID %>% str_split_fixed(., '_', 2))
  
  x$ID <- map(1:nrow(x), function(i) {
    y <- x[i, ]
    
    if(y[, 'X1']=='DQB1.0501' | y[, 'X1']=='DQB1.0602') {
      out <- y[, c('X1', 'X2')]
    } else out <- y[, c('X2', 'X1')]
    colnames(out) <- c('X1', 'X2')
    
    y <- data.frame(dplyr::select(y, -c(X1, X2)), out)
    
    y[, c('X1', 'X2')] %>% unlist %>% paste(collapse='_')
  }) %>% unlist
  
  x <- x[!(x$ID %>% duplicated), ] %>% 
    dplyr::select(-c(X1, X2))
  
  x$ID <- gsub('\\.', '*', x$ID) %>% 
    gsub('_', ' / ', .)
  
  return(data.frame(x, Pheno=p))
}

fread2 <- function(x) fread(x, data.table = F)


