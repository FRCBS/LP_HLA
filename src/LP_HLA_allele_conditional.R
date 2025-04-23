
## ===============================================
##
## Lichen planus HLA fine-mapping in R12 
## HLA allele conditional analyses
##
## ===============================================

source('src/LP_DQ_functions.R')


## -----------------------------------------------
## Covariate files for the cond. analyses
## -----------------------------------------------

makeNewCov(covars, 'DQB1*05:01') %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQB10501', sep='\t', quote=F, na='NA')

tmp <- fread('data/phenotypes/LP.covars.DQB10501', data.table = F)
tmp$DQB1.0501.domdev <- ifelse(tmp$DQB1.0501>0.5 & tmp$DQB1.0501<1.5, 1, 0)
fwrite(tmp, 'data/phenotypes/LP.covars.DQB10501.domdev', sep='\t', quote=F, na='NA') # 0501 and its domdev
fwrite(tmp %>% dplyr::select(-DQB1.0501), # 0501 domdev alone 
       'data/phenotypes/LP.covars.DQB10501.onlydomdev', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:01', 'DQB1*05:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10101_DQB10501', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10101', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:05')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10105', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:01', 'DQA1*01:05')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10101_DQA10105', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:05', 'DQB1*05:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10105_DQB10501', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:05', 'DQA1*01:01', 'DQB1*05:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10105_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:05', 'DQA1*01:01', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:05', 'DQA1*01:01', 'DQA1*01:02', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQA10105_DQA10101_DQA10102_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('B*35:01', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.B3501_DQA10105_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('B*13:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.B1302_DQA10105_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('B*08:01', 'B*13:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.B0801_B1302_DQA10105_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*03:01', 'B*08:01', 'B*13:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.A0301_B0801_B1302_DQA10105_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02', 'DRB1*09:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02', 'DRB1*09:01', 'DQB1*02:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02', 'DRB1*09:01', 'DQB1*02:02', 'B*07:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*03:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.A0301', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQB1*06:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQB10602', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'B*07:02')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.Class1', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02', 'DRB1*09:01', 'DQB1*02:02', 'DQB1*03:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.Class2', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('A*03:01', 'A*02:01', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02', 'DRB1*09:01', 'DQB1*02:02', 'DQB1*03:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.A_Class2', sep='\t', quote=F, na='NA')

makeNewCov(covars, c('B*08:01', 'B*13:02', 'B*07:02', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*06:02', 'DRB1*09:01', 'DQB1*02:02', 'DQB1*03:01')) %>% 
  fwrite(., 'data/phenotypes/LP.covars.B_Class2', sep='\t', quote=F, na='NA')

tmp <- makeNewCov(covars, c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'B*07:02', 'DQA1*01:05', 'DQA1*01:01', 'DQB1*05:01'))
tmp$DQB10501hom <- 0
tmp$DQB10501hom[(tmp[, 'DQA1.0105'] > 1.5 & tmp[, 'DQB1.0501'] > 1.5) | 
                  (tmp[, 'DQA1.0101'] > 1.5 & tmp[, 'DQB1.0501'] > 1.5)] <- 1
tmp %>% dplyr::select(-c(DQA1.0105, DQA1.0101, DQB1.0501)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.Class1_DQB10501hom', sep='\t', quote=F, na='NA')
tmp %>% dplyr::select(-c(A.0201, A.0301, B.0801, B.1302, B.0702, DQA1.0105, DQA1.0101, DQB1.0501)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQB10501hom', sep='\t', quote=F, na='NA')
tmp$DQB10501het <- 0
tmp$DQB10501het[(tmp[, 'DQA1.0105'] > 0.5 & tmp[, 'DQB1.0501'] > 0.5 & tmp$DQB10501hom==0) | 
                  (tmp[, 'DQA1.0101'] > 0.5 & tmp[, 'DQB1.0501'] > 0.5  & tmp$DQB10501hom==0)] <- 1
tmp %>% dplyr::select(-c(DQA1.0105, DQA1.0101, DQB1.0501, DQB10501hom)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.Class1_DQB10501het', sep='\t', quote=F, na='NA')
tmp %>% dplyr::select(-c(A.0201, A.0301, B.0801, B.1302, B.0702, DQA1.0105, DQA1.0101, DQB1.0501, DQB10501hom)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQB10501het', sep='\t', quote=F, na='NA')

tmp <- makeNewCov(covars, c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'B*07:02', 'DQA1*01:02', 'DQB1*06:02'))
tmp$DQB10602hom <- 0
tmp$DQB10602hom[(tmp[, 'DQA1.0102'] > 1.5 & tmp[, 'DQB1.0602'] > 1.5)] <- 1
tmp %>% dplyr::select(-c(DQA1.0102, DQB1.0602)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.Class1_DQB10602hom', sep='\t', quote=F, na='NA')
tmp %>% dplyr::select(-c(A.0201, A.0301, B.0801, B.1302, B.0702, DQA1.0102, DQB1.0602)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQB10602hom', sep='\t', quote=F, na='NA')
tmp$DQB10602het <- 0
tmp$DQB10602het[(tmp[, 'DQA1.0102'] > 0.5 & tmp[, 'DQB1.0602'] > 0.5 & tmp$DQB10602hom==0)] <- 1
tmp %>% dplyr::select(-c(DQA1.0102, DQB1.0602, DQB10602hom)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.Class1_DQB10602het', sep='\t', quote=F, na='NA')
tmp %>% dplyr::select(-c(A.0201, A.0301, B.0801, B.1302, B.0702, DQA1.0102, DQB1.0602, DQB10602hom)) %>% 
  fwrite(., 'data/phenotypes/LP.covars.DQB10602het', sep='\t', quote=F, na='NA')





## -----------------------------------------------
## Genotype file for 0501 domdev
## -----------------------------------------------

tmp <- fread('data/phenotypes/LP.covars.DQB10501', data.table = F)
tmp$DQB1.0501.domdev <- ifelse(tmp$DQB1.0501>0.5 & tmp$DQB1.0501<1.5, 1, 0)

createPlinkGeno(tmp %>% dplyr::select(c(IID, DQB1.0501.domdev)), '0501domdev')



## -----------------------------------------------
## regenie step 1: null models
## v3.0.1
## -----------------------------------------------

# Null model with basic covariates and all 3 phenotypes
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              " --extract data/LP_regenie_null_qc_pass.snplist ", 
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars ",
              " --bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              " --out results/LP_DQ/hm3_step1_Basic_All"))

# DQB1*0501 dosage as an additional covariate
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10501 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10501"))

# DQB1*0501 dosage and domdev as an additional covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10501.domdev ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10501_domdev"))

# DQB1*0501 domdev alone as an additional covariate
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10501.onlydomdev ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10501_onlydomdev"))

# DQA1*0101 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10101 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10101"))

# DQA1*0105 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10105 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10105"))

# DQA1*0101, DQA1*0105 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10101_DQA10105 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10101_DQA10105"))

# DQA1*0101 and DQB1*0501 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10101_DQB10501 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10101_DQB10501"))

# DQA1*0105 and DQB1*0501 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10105_DQB10501 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10105_DQB10501"))

# DQA1*0105, DQA1*0101, and DQB1*0501 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQB10501"))

# DQA1*0105 and DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10105_DQB10501_DQB10602"))

# DQA1*0105, DQA1*0101, DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQB10501_DQB10602"))

# DQA1*0105, DQA1*0101, DQA1*0102, DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQA10102_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQA10102_DQB10501_DQB10602"))

# B*3501, DQA1*0105 and DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.B3501_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_B3501_DQA10105_DQB10501_DQB10602"))

# B*1302, DQA1*0105 and DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_B1302_DQA10105_DQB10501_DQB10602"))

# B*0801, B*1302, DQA1*0105 and DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.B0801_B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_B0801_B1302_DQA10105_DQB10501_DQB10602"))

# A*0301, B*0801, B*1302, DQA1*0105 and DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.A0301_B0801_B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602"))

# A*0201, A*0301, B*0801, B*1302, DQA1*0105 and DQB1*0501 and DQB1*0602 as additional additive covariates
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602"))

# A*0201, A*0301, B*0801, B*1302, DQA1*0105 and DQB1*0501, DQB1*0602 and DRB1*0901 as additional additive covariates
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901"))

# A*0201, A*0301, B*0801, B*1302, DQA1*0105 and DQB1*0501, DQB1*0602, DRB1*0901 and DQB1*0202 as additional additive covariates
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202"))

# A*0201, A*0301, B*0801, B*1302, B*0702, DQA1*0105 and DQB1*0501, DQB1*0602, DRB1*0901 and DQB1*0202 as additional additive covariates
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202"))

# A*0301
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.A0301 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_A0301"))

# DQB1*0602
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10602 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10602"))

# Class1
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_Class1"))

# Class2
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class2 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_Class2"))

# A and Class2
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.A_Class2 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_A_Class2"))

# B and Class2
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.B_Class2 ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_B_Class2"))

# Class I and 0501 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10501hom ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_Class1_DQB10501hom"))

# 0501 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10501hom ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10501hom"))

# Class I and 0602 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10602hom ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_Class1_DQB10602hom"))

# 0602 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10602hom ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10602hom"))

# Class I and 0501 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10501het ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_Class1_DQB10501het"))

# 0501 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10501het ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10501het"))

# Class I and 0602 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10602het ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_Class1_DQB10602het"))

# 0602 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              "--extract data/LP_regenie_null_qc_pass.snplist --keep data/LP_regenie_null_qc_pass.id ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10602het ",
              "--bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              "--out results/LP_DQ/hm3_step1_DQB10602het"))



## -----------------------------------------------
## regenie step 2: association
## v3.0.1
## -----------------------------------------------

# basic = HLA association without any additional allele covariates
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Basic_All_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_basic"))

# HLA association adjusted for DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10501"))

# HLA association adjusted for DQB1*0501 dosages and domdev
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQB10501_domdev_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10501_domdev"))

# 0501 domdev association adjusted for DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bed data/genotypes/0501domdev_haplos ",
              "--pred results/LP_DQ/hm3_step1_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_0501domdev_step2_DQB10501"))

# HLA association adjusted for DQB1*0501 domdev alone
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQB10501_onlydomdev_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10501.onlydomdev ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10501_onlydomdev"))

# HLA association adjusted for DQA1*0101 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10101_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10101 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10101"))

# HLA association adjusted for DQA1*0105 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10105"))

# HLA association adjusted for DQA1*0101, DQA1*0105 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10101_DQA10105_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10101_DQA10105 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10101_DQA10105"))

# HLA association adjusted for DQA1*0101, DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10101_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10101_DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10101_DQB10501"))

# HLA association adjusted for DQA1*0105, DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10105_DQB10501"))

# HLA association adjusted for DQA1*0105, DQA1*0101, DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10105_DQA10101_DQB10501"))

# HLA association adjusted for DQA1*0105, DQA1*0101, DQB1*0501, DQB1*0602 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10105_DQA10101_DQB10501_DQB10602"))

# HLA association adjusted for DQA1*0105, DQA1*0101, DQA1*0102, DQB1*0501, DQB1*0602 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQA10102_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQA10102_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10105_DQA10101_DQA10102_DQB10501_DQB10602"))

# association adjusted for DQA1*0105, DQB1*0602 and DQB1*0501
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQA10105_DQB10501_DQB10602"))

# association adjusted for B*3501, DQA1*0101, DQB1*0602 and DQB1*0501 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_B3501_DQA10105_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.B3501_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_B3501_DQA10105_DQB10501_DQB10602"))

# association adjusted for B*1302, DQA1*0105, DQB1*0602 and DQB1*0501 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_B1302_DQA10105_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_B1302_DQA10105_DQB10501_DQB10602"))

# association adjusted for B*0801, B*1302, DQA1*0105, DQB1*0602 and DQB1*0501 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_B0801_B1302_DQA10105_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.B0801_B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_B0801_B1302_DQA10105_DQB10501_DQB10602"))

# association adjusted for A*0301, B*0801, B*1302, DQA1*0105, DQB1*0602 and DQB1*0501 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0301_B0801_B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602"))

# association adjusted for A*0201, A*0301, B*0801, B*1302, DQA1*0105, DQB1*0602 and DQB1*0501 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602"))

# association adjusted for A*0201, A*0301, B*0801, B*1302, DQA1*0105, DQB1*0602, DQB1*0501 and DRB1*0901 
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901"))

# association adjusted for A*0201, A*0301, B*0801, B*1302, DQA1*0105, DQB1*0602, DQB1*0501, DRB1*0901, DQB1*0202 
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202"))

# association adjusted for A*0201, A*0301, B*0801, B*1302, B*0702, DQA1*0105, DQB1*0602, DQB1*0501, DRB1*0901, DQB1*0202 
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202"))

# A*0301
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0301_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0301 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_A0301"))

# DQB1*0602
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10602"))


# Class I 0501 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class1_DQB10501hom_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.Class1_DQB10501hom ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_Class1_DQB10501hom"))

# 0501 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQB10501hom_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10501hom ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10501hom"))

# Class I 0602 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class1_DQB10602hom_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.Class1_DQB10602hom ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_Class1_DQB10602hom"))

# 0602 hom
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class1_DQB10602hom_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10602hom ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10602hom"))

# Class I 0501 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class1_DQB10501het_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.Class1_DQB10501het ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_Class1_DQB10501het"))

# 0501 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQB10501het_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10501het ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10501het"))

# Class I 0602 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class1_DQB10602het_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.Class1_DQB10602het ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_Class1_DQB10602het"))

# 0602 het
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class1_DQB10602het_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10602het ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_step2_DQB10602het"))



## -----------------------------------------------
## results
## -----------------------------------------------

# 0501 domdev assoc
rbind(
  process_result('results/LP_DQ/LP_0501domdev_step2_DQB10501_All.regenie', 'All'),
  process_result('results/LP_DQ/LP_0501domdev_step2_DQB10501_NonOral.regenie', 'NonOral'),
  process_result('results/LP_DQ/LP_0501domdev_step2_DQB10501_Oral.regenie', 'Oral')
) %>% dplyr::select(c(ID, A1FREQ, BETA, SE, P, Pheno)) %>% kable(format='simple') 

# HLA adj for 0501 and domdev
rbind(
  process_result('results/LP_DQ/LP_HLA_step2_DQB10501_domdev_All.regenie', 'All'),
  process_result('results/LP_DQ/LP_HLA_step2_DQB10501_domdev_NonOral.regenie', 'NonOral'),
  process_result('results/LP_DQ/LP_HLA_step2_DQB10501_domdev_Oral.regenie', 'Oral')
) %>% dplyr::select(c(ID, A1FREQ, BETA, SE, P, Pheno)) %>% filter(ID == 'DQB1*05:01') %>% 
  kable(format='simple') 

# HLA adj for 0501domdev alone
rbind(
  process_result('results/LP_DQ/LP_HLA_step2_DQB10501_onlydomdev_All.regenie', 'All'),
  process_result('results/LP_DQ/LP_HLA_step2_DQB10501_onlydomdev_NonOral.regenie', 'NonOral'),
  process_result('results/LP_DQ/LP_HLA_step2_DQB10501_onlydomdev_Oral.regenie', 'Oral')
) %>% dplyr::select(c(ID, A1FREQ, BETA, SE, P, Pheno)) %>% filter(ID == 'DQB1*05:01') %>% 
  kable(format='simple') 

# just HLA
rbind(
  process_result('results/LP_DQ/LP_HLA_step2_basic_All.regenie', 'All'),
  process_result('results/LP_DQ/LP_HLA_step2_basic_NonOral.regenie', 'NonOral'),
  process_result('results/LP_DQ/LP_HLA_step2_basic_Oral.regenie', 'Oral')
) %>% dplyr::select(c(ID, A1FREQ, BETA, SE, P, Pheno)) %>% filter(ID == 'DQB1*05:01') %>% 
  kable(format='simple') 

# in temporal order
system("ls -t results/LP_DQ/LP_HLA_step2_*_All.regenie")

# list regenie assoc data for the adjustment chain
res <- list(list.files('results/LP_DQ', 'LP_HLA_step2_basic_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQB10501_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQB10501_domdev_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQA10101_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQA10105_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQA10101_DQA10105_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQA10101_DQB10501_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQA10105_DQB10501_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_DQA10105_DQB10501_DQB10602_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_B3501_DQA10105_DQB10501_DQB10602_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_B1302_DQA10105_DQB10501_DQB10602_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_B0801_B1302_DQA10105_DQB10501_DQB10602_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_step2_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_[A|O|N]+', full.names=T))

# adjustment alleles in the same order
adjustment.chain <- c('none', 
                      'DQB1*0501', 
                      'DQB1*0501 domdev', 
                      'DQA1*0101', 
                      'DQA1*0105', 
                      'DQA1*0101 & DQA1*0105', 
                      'DQA1*0101 & DQB1*0501', 
                      'DQA1*0105 & DQB1*0501', 
                      'DQA1*0105 & DQB1*0501 &\nDQB1*0602',
                      'DQA1*0105 & DQB1*0501 &\nDQB1*0602 & B*3501', 
                      'DQA1*0105 & DQB1*0501 &\nDQB1*0602 & B*1302',
                      'DQA1*0105 & DQB1*0501 &\nDQB1*0602 & B*1302 & B*0801',
                      'DQA1*0105 & DQB1*0501 &\nDQB1*0602 & B*1302 & B*0801 &\nA*0201 & A*0301 & DRB1*0901')

# interesting alleles 'DQA1*02:01', 'DQB1*02:01',  'DQA1*01:01',  'DQA1*01:05', 
target.alleles <- c('A*03:01', 'B*13:02', 'B*35:01', 'B*07:02', 'B*08:01', 
                    'DQA1*01:02', 
                    'DQB1*03:01', 'DQB1*02:01', 'DQB1*05:01', 'DQB1*02:02', 'DQB1*06:02', 
                    'DRB1*09:01')
# DQB1*0501 with and without domdev
# All
plot(fread(res[[2]][[1]])$BETA, fread(res[[3]][[1]])$BETA, 
     xlab = 'DQB1*0501', ylab = 'DQB1*0501 domdev', main = 'All')
plot(fread(res[[2]][[2]])$BETA, fread(res[[3]][[2]])$BETA, 
     xlab = 'DQB1*0501', ylab = 'DQB1*0501 domdev', main = 'NonOral')
plot(fread(res[[2]][[3]])$BETA, fread(res[[3]][[3]])$BETA, 
     xlab = 'DQB1*0501', ylab = 'DQB1*0501 domdev', main = 'Oral')

# read allele results data
res.alleles <- map_dfr(1:length(res), function(i) {
  out <- map2(res[[i]], c('All', 'NonOral', 'Oral'), function(d, p) process_result(d, p)) %>% 
    do.call(rbind, .) %>% filter(ID %in% target.alleles) %>% 
    arrange(P)
  out$ID <- gsub(':', '', out$ID, fixed=T)
  out$Adjusted <- adjustment.chain[i]
  out$Signif <- ifelse(out$FDR_BY<0.05, 'FDR <0.05', 'not significant')
  out %>% return()
})
# select adjustments
res.alleles <- res.alleles %>% filter(Adjusted %in% adjustment.chain)#[c(1, 3, 4, 6)])
res.alleles$Adjusted <- factor(res.alleles$Adjusted, levels=adjustment.chain)#[c(1, 3, 4, 6)])
# arrange alleles by beta
res.alleles$ID <- factor(res.alleles$ID, levels=res.alleles %>% filter(Adjusted=='none') %>%  
                           group_by(ID) %>% summarise(MB=mean(BETA)) %>% arrange(MB) %>% .$ID)
# arrange phenos by beta
res.alleles$Pheno <- factor(res.alleles$Pheno, levels=c('Oral', 'All', 'NonOral'))

# limit to  DQA1, DQB1
target.alleles.2 <- c('DQA1*01:01', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*03:01',
                      'B*08:01', 'B*13:02', 'B*35:01', 'DQA1*01:02', 'DQB1*06:02')

# read allele results data
res.alleles.2 <- map_dfr(1:length(res), function(i) {
  out <- map2(res[[i]], c('All', 'NonOral', 'Oral'), function(d, p) process_result(d, p)) %>% 
    do.call(rbind, .) %>% filter(ID %in% target.alleles.2) %>% 
    arrange(P)
  out$ID <- gsub(':', '', out$ID, fixed=T)
  out$Adjusted <- adjustment.chain[i]
  out$Signif <- ifelse(out$FDR_BY<0.05, 'FDR <0.05', 'not significant')
  out %>% return()
})
# select adjustments
res.alleles.2 <- res.alleles.2 %>% filter(Adjusted %in% adjustment.chain[c(1:2, 6:7)])
res.alleles.2$Adjusted <- factor(res.alleles.2$Adjusted, levels=adjustment.chain[c(1:2, 6:7)])
# arrange alleles by beta
res.alleles.2$ID <- factor(res.alleles.2$ID, levels=res.alleles.2 %>% filter(Adjusted=='none') %>%  
                           group_by(ID) %>% summarise(MB=mean(BETA)) %>% arrange(MB) %>% .$ID)
# arrange phenos by beta
res.alleles.2$Pheno <- factor(res.alleles.2$Pheno, levels=c('Oral', 'All', 'NonOral'))


# limit to  DQA1 adjustments
target.alleles.3 <- c('DQA1*01:01', 'DQA1*01:05', 'DQB1*05:01', 'DQB1*03:01',
                      'B*08:01', 'B*13:02', 'B*35:01', 'DQA1*01:02', 'DQB1*06:02')

# read allele results data
res.alleles.3 <- map_dfr(1:length(res), function(i) {
  out <- map2(res[[i]], c('All', 'NonOral', 'Oral'), function(d, p) process_result(d, p)) %>% 
    do.call(rbind, .) %>% filter(ID %in% target.alleles.3) %>% 
    arrange(P)
  out$ID <- gsub(':', '', out$ID, fixed=T)
  out$Adjusted <- adjustment.chain[i]
  out$Signif <- ifelse(out$FDR_BY<0.05, 'FDR <0.05', 'not significant')
  out %>% return()
})
# select adjustments
res.alleles.3 <- res.alleles.3 %>% filter(Adjusted %in% adjustment.chain[c(3:5)])
res.alleles.3$Adjusted <- factor(res.alleles.3$Adjusted, levels=adjustment.chain[c(3:5)])
# arrange alleles by beta
res.alleles.3$ID <- factor(res.alleles.3$ID, levels=res.alleles.3  %>%  
                             group_by(ID) %>% summarise(MB=mean(BETA)) %>% arrange(MB) %>% .$ID)
# arrange phenos by beta
res.alleles.3$Pheno <- factor(res.alleles.3$Pheno, levels=c('Oral', 'All', 'NonOral'))


## -----------------------------------------------
## allele correlation plot
## -----------------------------------------------

jpeg('results/LP_DQ/Allele_correlation.jpg', width = 7, height = 7, res = 1000, units = 'in')
hla %>% dplyr::select(all_of(
  c('A*02:01', 'A*03:01', 'B*08:01', 'B*13:02', 'B*35:01', 'B*07:02',  
    'DQA1*01:01', 'DQA1*01:02', 'DQA1*01:05',
    'DQB1*02:01', 'DQB1*02:02', 'DQB1*03:01', 'DQB1*05:01', 'DQB1*06:02', 
    'DRB1*09:01'))) %>% rename_with(., function(x) gsub(':', '', x)) %>% 
  cor %>% corrplot(tl.col = 'grey20')
dev.off()


## -----------------------------------------------
## print top tables
## -----------------------------------------------

# table 1 of summary stats
i <- 1
out <- map2(res[[i]], c('All', 'NonOral', 'Oral'), function(d, p) process_result(d, p)) %>% 
  do.call(rbind, .) %>% 
  arrange(P)
out$ID <- gsub(':', '', out$ID, fixed=T)
out.ids <- out %>% filter(Pheno == 'All') %>% head(30) %>% .$ID
out %>% filter(ID %in% out.ids) %>% dplyr::select(c(ID, A1FREQ, Pheno, P, OR, OR_95CI)) %>% 
  pivot_wider(names_from = Pheno, values_from = c(P, OR, OR_95CI))

list(
  out %>% filter(ID %in% out.ids, Pheno == 'All') %>% 
    dplyr::select(c(ID, A1FREQ, Pheno, P, OR, OR_95CI)),
  out %>% filter(ID %in% out.ids, Pheno == 'NonOral') %>% 
    dplyr::select(c(ID, A1FREQ, Pheno, P, OR, OR_95CI)),
  out %>% filter(ID %in% out.ids, Pheno == 'Oral') %>% 
    dplyr::select(c(ID, A1FREQ, Pheno, P, OR, OR_95CI))) %>% 
  Reduce(function(x, y) inner_join(x, y, by = 'ID'), .) %>% 
  fwrite(., 'results/LP_DQ/Table1.tsv', sep = '\t', quote = F)

# table 2 of allele freqs
computeFreq <- function(x) apply(x, 2, function(k) sqrt(sum(k>1.5)/length(k))) %>% data.frame

freqs.all <- inner_join(
  hla %>% filter(IID %in% filter(phenos, All==1)$IID) %>% dplyr::select(all_of(target.alleles)) %>% 
    computeFreq %>% rownames_to_column(),
  hla %>% filter(IID %in% filter(phenos, All==0)$IID) %>% dplyr::select(all_of(target.alleles)) %>% 
    computeFreq %>% rownames_to_column(), by = 'rowname')
colnames(freqs.all) <- c('HLA', 'All', 'Controls') 

freqs.oral <- inner_join(
  hla %>% filter(IID %in% filter(phenos, Oral==1)$IID) %>% dplyr::select(all_of(target.alleles)) %>% 
    computeFreq %>% rownames_to_column(),
  hla %>% filter(IID %in% filter(phenos, Oral==0)$IID) %>% dplyr::select(all_of(target.alleles)) %>% 
    computeFreq %>% rownames_to_column(), by = 'rowname')
colnames(freqs.oral) <- c('HLA', 'Oral', 'Controls') 

freqs.nonoral <- inner_join(
  hla %>% filter(IID %in% filter(phenos, NonOral==1)$IID) %>% dplyr::select(all_of(target.alleles)) %>% 
    computeFreq %>% rownames_to_column(),
  hla %>% filter(IID %in% filter(phenos, NonOral==0)$IID) %>% dplyr::select(all_of(target.alleles)) %>% 
    computeFreq %>% rownames_to_column(), by = 'rowname')
colnames(freqs.nonoral) <- c('HLA', 'NonOral', 'Controls') 

Reduce(function(x, y) inner_join(x, y, by = 'HLA'), 
       list(freqs.all, freqs.nonoral, freqs.oral)) %>% 
  fwrite('results/LP_DQ/Table2_Allele_frequencies.tsv', sep = '\t', quote = F)


# addition to table 2: HLA haplotype freqs

target.haplotypes <- list(c('DQA1*01:01', 'DQB1*05:01'),
                          c('DQA1*01:05', 'DQB1*05:01'),
                          c('DRB1*15:01', 'DQA1*01:02', 'DQB1*06:02'))
hap <- map(target.haplotypes, function(x) {
  # x <- target.haplotypes[[1]]
  tmp <- dplyr::select(hla, all_of(x))
  if(ncol(tmp)==3) {
    tmp.ind <- which(tmp[, 1]>1.5 & tmp[, 2]>1.5 & tmp[, 3]>1.5)
  } else tmp.ind <- which(tmp[, 1]>1.5 & tmp[, 2]>1.5)
  tmp.res <- data.frame(IID = hla$IID, HAP = 0)
  tmp.res$HAP[tmp.ind] <- 1
  colnames(tmp.res)[2] <- paste(x, collapse='_') %>% 
    gsub('\\*', '.', .) %>% gsub('\\:', '', .)
  tmp.res %>% return()
}) %>% Reduce(inner_join, .)

computeFreq <- function(x) apply(x, 2, function(k) sqrt(sum(k==1)/length(k))) %>% data.frame

freqs.all <- inner_join(
  hap %>% filter(IID %in% filter(phenos, All==1)$IID) %>% dplyr::select(-IID) %>% 
    computeFreq %>% rownames_to_column(),
  hap %>% filter(IID %in% filter(phenos, All==0)$IID) %>% dplyr::select(-IID) %>%  
    computeFreq %>% rownames_to_column(), by = 'rowname')
colnames(freqs.all) <- c('hap', 'All', 'Controls') 

freqs.oral <- inner_join(
  hap %>% filter(IID %in% filter(phenos, Oral==1)$IID) %>%  dplyr::select(-IID) %>% 
    computeFreq %>% rownames_to_column(),
  hap %>% filter(IID %in% filter(phenos, Oral==0)$IID) %>%  dplyr::select(-IID) %>% 
    computeFreq %>% rownames_to_column(), by = 'rowname')
colnames(freqs.oral) <- c('hap', 'Oral', 'Controls') 

freqs.nonoral <- inner_join(
  hap %>% filter(IID %in% filter(phenos, NonOral==1)$IID) %>%  dplyr::select(-IID) %>% 
    computeFreq %>% rownames_to_column(),
  hap %>% filter(IID %in% filter(phenos, NonOral==0)$IID) %>%  dplyr::select(-IID) %>% 
    computeFreq %>% rownames_to_column(), by = 'rowname')
colnames(freqs.nonoral) <- c('hap', 'NonOral', 'Controls') 

Reduce(function(x, y) inner_join(x, y, by = 'hap'), 
       list(freqs.all, freqs.nonoral, freqs.oral)) %>% 
  fwrite('results/LP_DQ/Table2.1_Haplo_frequencies.tsv', sep = '\t', quote = F)


  
## -----------------------------------------------
## plots
## -----------------------------------------------

# filter results for plotting
tmp <- res.alleles %>% 
  filter(!(Adjusted %in% c('DQA1*0101', 
                           'DQA1*0105', 
                           'DQA1*0101 & DQB1*0501',
                           'DQA1*0101 & DQA1*0105',
                           'DQB1*0501 domdev'))) %>%
  filter(Pheno %in% c('All', 'Oral', 'NonOral')) 
tmp$Adjusted %>% as.character %>% unique

# panel annotation
panel.anno <- data.frame(OR = 0.7, ID = 13.8, Pheno = tmp$Pheno[1], OR_LOWER = 1, OR_UPPER = 1, 
                         Signif = tmp$Signif[1],
                         label = rep(paste0(letters[1:8], ')')),
                         Adjusted = tmp$Adjusted %>% as.character %>% unique %>% factor) 

p.res.alleles <- ggplot(tmp, 
                        aes(OR, ID, xmin = OR_LOWER, xmax = OR_UPPER, color = Pheno %>% factor, shape = Signif)) +
  geom_pointrange(position=position_dodge(width=0.6), size=0.6, linewidth=0.7) +
  xlab('OR') + ylab('') +
  scale_color_manual(values=c("#440154FF", "#21908CFF",  "#FFBF00") %>% rev) + # , #FFD580 ff9248
  scale_shape_manual(name=NULL, values=c(15, 0)) +
  scale_x_continuous(transform = 'log10') +
  geom_text(data = panel.anno, aes(OR, ID, label = label, xmin = OR_LOWER, xmax = OR_UPPER, shape = Signif), 
            size = 4.3, color = 'black', fontface = 'bold') +
  coord_cartesian(xlim=c(0.7, 2), ylim = c(1, 12), clip = 'off') +
  geom_vline(xintercept=1, linewidth=0.3, color='black', linetype='dashed') +
  guides(color = guide_legend(override.aes=list(shape = 15))) +
  facet_wrap(vars(Adjusted), scales='fixed', ncol=3, nrow = 3) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, angle=0, hjust=0.5),
        axis.text.y=element_text(size=10),
        panel.grid=element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.line=element_blank(),
        legend.position='inside',
        legend.position.inside = c(0.82, 0.13),
        legend.title=element_blank(),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
        panel.border=element_rect(fill=NA, linewidth=0.2))
sapply(seq(1, unique(res.alleles$ID) %>% length, by=2), function(i) {
  p.res.alleles <<- p.res.alleles +
    annotate(geom='rect', ymin=i-0.5, ymax=i+0.5, xmin= 0, xmax=Inf, 
             fill='grey25', alpha=0.05, color=NA) 
}) 

# Figure 3
jpeg('results/LP_DQ/Fig3_Alleles_adj.jpg', width=10, height=11, units='in', res=1000)
p.res.alleles
dev.off()
