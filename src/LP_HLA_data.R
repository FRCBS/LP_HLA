
## ===============================================
##
## Lichen planus HLA fine-mapping in R12 
## Data read-in and processing
##
## ===============================================

source('src/LP_DQ_functions.R')

## -----------------------------------------------
## HLA dosages
## -----------------------------------------------

# # HLA BGEN file
# hla.bgen   <- '/finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen'
# hla.sample <- '/finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample'
# 
# # convert to plink dosage format
# system(paste0("plink2 --bgen ", hla.bgen, " ref-first --sample ", hla.sample,
#               " --sort-vars --make-pgen --out results/dosages/R12_HLA"))
# system(paste0("plink2 --pfile results/dosages/R12_HLA --recode A --out results/dosages/R12_HLA"))

# read HLA dosages
hla <- fread('results/dosages/R12_HLA.raw', data.table=F)[, -c(3:6)]
hla[, -c(1:2)] <- 2 - hla[, -c(1:2)] # 
colnames(hla) <- gsub('_<absent>', '', colnames(hla), fixed=T)

## -----------------------------------------------
## Covariate data
## -----------------------------------------------

covars <- fread('/finngen/library-red/finngen_R12/analysis_covariates/R12_COV_V2.FID.txt.gz', data.table=F)
covars$regionofbirthname <- gsub('-| ', '', covars$regionofbirthname) %>% 
  gsub('ä', 'a', .) %>% gsub('Å', 'A', .)
batches <- colnames(covars)[covars %>% colnames %>% grepl('BATCH_DS', .)]
covars <- dplyr::select(covars, c(FID, IID, AGE_AT_DEATH_OR_END_OF_FOLLOWUP, SEX, regionofbirthname,
                                  IS_FINNGEN2_CHIP, one_of(batches), BL_YEAR,
                                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))
covars <- covars %>% na.omit
covars %>% dim() # n=499483
covars %>% head

# categorical variables to numeric
covars <- data.frame(covars[, 1:2], model.matrix( ~ ., data=covars[, -c(1:2)]))
covars <- dplyr::select(covars, -X.Intercept.)
covars %>% dim # n=499483      
fwrite(covars, 'data/phenotypes/LP.covars', sep='\t', quote=F)
covars <- fread('data/phenotypes/LP.covars')


## -----------------------------------------------
## Phenotype data
## -----------------------------------------------

# # subtype sample ID lists
# system(paste0('tar -xvf /finngen/shared/lp_subtype_lists/20230605_230001/files/mpreeve/LP_R11.tar -C ', 
#               'data/phenotypes'))

lp.subtypes <- c('AllLP.tsv', 'OralLP.tsv', 'NonOralLP.tsv') %>% 
  paste0('data/phenotypes/', .) %>% 
  map(fread, data.table=F, header=T) %>% Reduce(full_join, .)
lp.subtypes %>% dim()
colnames(lp.subtypes)[-1] <- c('All', 'Oral', 'NonOral')
lp.subtypes %>% filter(All==1) %>% dim # 8484
lp.subtypes %>% filter(Oral==1) %>% dim # 3663
lp.subtypes %>% filter(NonOral==1) %>% dim # 4821

# write pheno file for Regenie
fwrite(data.frame(FID=lp.subtypes$IID, lp.subtypes), 'data/phenotypes/LP.pheno', sep='\t', na='NA', quote=F)

# read processed pheno data
phenos <- fread('data/phenotypes/LP.pheno')
phenos %>% head
