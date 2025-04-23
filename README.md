# LP_HLA

### Analysis scripts for HLA fine-mapping for oral and non-oral lichen planus in FinnGen. 

J Ritari, MP Reeve, FinnGen, M Siponen, M Vehviläinen, T Salo, K Osoegawa, M Fernandez Viña, B Goudey, J Partanen, E Mignot. Fine mapping of HLA effects in Oral and non-Oral lichen planus. _Submitted manuscript._

Analysis scripts in `src/`:

`LP_HLA_allele_combinations.R`
Analysis of known combinations of HLA-DQA1 and HLA-DQB1 alleles as they occur in FinnGen.
Definition of DQ allele pairs and their dosages. LP association analysis using Regenie.
Plotting the results.

`LP_HLA_allele_conditional.R`
Conditional association analysis of HLA alleles using Regenie. Specific HLA alleles are included as additional model covariates. Plotting the results.

`LP_HLA_data.R`
Data preparation steps prior to analysis. Conversion of HLA genotype data to dosages. Preparation of basic covariate and LP phenotype matrices for use with Regenie.

`LP_HLA_functions.R`
Helper functions for analyses. Loading R libraries, functions for converting dosage information to plink format for input to Regenie. Functions for processing association output from Regenie.

`LP_HLA_homozygotes.R`
Definition and extraction of DQA1-DQB1 homozygotes from data. HLA allele association analayses on homozygotic samples using Regenie. Plotting the results.

`LP_HLA_SNPs.R`
Association analyses on MHC region SNP genotypes conditioned on specific HLA alleles using Regenie. Plotting the results.
