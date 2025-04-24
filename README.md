# LP_HLA

### Analysis scripts for HLA fine-mapping for oral and non-oral lichen planus in FinnGen. 


J Ritari, MP Reeve, FinnGen, M Siponen, M Vehviläinen, T Salo, K Osoegawa, M Fernandez Viña, B Goudey, J Partanen, E Mignot. **Fine mapping of HLA effects in Oral and non-Oral lichen planus.** _Submitted manuscript._


#### Scripts in `src/`:

+ `LP_HLA_allele_combinations.R`
Analysis of known combinations of HLA-DQA1 and HLA-DQB1 alleles as they occur in imputed FinnGen HLA data. Definition of DQ allele pairs and their dosages. DQ combination genotype association analysis using Regenie. Plotting the results.

+ `LP_HLA_allele_conditional.R`
Conditional association analysis of HLA alleles using Regenie. Specific HLA allele dosages are included as additional model covariates. Plotting the results.

+ `LP_HLA_data.R`
Data preparation steps prior to analysis. Conversion of imputed HLA genotype data to dosages. Preparation of basic covariate and LP phenotype matrices for use with Regenie.

+ `LP_HLA_functions.R`
Helper functions for analyses. Loading R libraries, functions for converting dosage information to plink format for input to Regenie. Functions for processing association output files produced by Regenie.

+ `LP_HLA_homozygotes.R`
Definition and extraction of DQA1-DQB1 homozygotes from imputed HLA allele data. HLA allele association analayses on DQ homozygotic samples using Regenie. Plotting the results.

+ `LP_HLA_SNPs.R`
Association analyses on MHC region SNP genotypes conditioned on specific HLA allele dosages using Regenie. Plotting the results.
