
## ===============================================
##
## Lichen planus HLA fine-mapping in R12 
## HLA allele conditional analyses
##
## ===============================================

source('src/LP_DQ_functions.R')


## -----------------------------------------------
## Extract MHC SNPs
## -----------------------------------------------

system(paste0("plink2 --bfile /finngen/library-red/finngen_R12/genotype_plink_3.0/data/finngen_R12 ", 
              "--chr 6 --from-mb 29 --to-mb 34 --make-pgen --out ./data/genotypes/R12_MHC"))


## -----------------------------------------------
## plink assoc for MHC SNPs
## condition with class II, class I, both
## -----------------------------------------------

system(paste0("plink2 --pfile data/genotypes/R12_MHC --pheno data/phenotypes/LP.pheno --1 ", 
              "--covar data/phenotypes/LP.covars --glm hide-covar --covar-variance-standardize ",
              "--out results/LP_DQ/plink/basic"))
              



## -----------------------------------------------
## Full HLA adjustment for MHC SNPs
## -----------------------------------------------


# run regenie
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --pgen data/genotypes/R12_MHC ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_SNP_step2_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --pgen data/genotypes/R12_MHC ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_SNP_step2_A0201_A0301_B0801_B1302_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --pgen data/genotypes/R12_MHC ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_SNP_step2_A0201_A0301_B0801_B1302_B0702_DQA10105_DQB10501_DQB10602_DRB10901_DQB10202"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --pgen data/genotypes/R12_MHC ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class1_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.Class1 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_SNP_step2_Class1"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --pgen data/genotypes/R12_MHC ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Class2_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.Class2 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_SNP_step2_Class2"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --pgen data/genotypes/R12_MHC ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_A_Class2_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.A_Class2 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_SNP_step2_A_Class2"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              "--step 2 --pgen data/genotypes/R12_MHC ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_B_Class2_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.B_Class2 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_SNP_step2_B_Class2"))


## -----------------------------------------------
## Results
## -----------------------------------------------

ls.snp <- list.files('results/LP_DQ', 'SNP_step2.*regenie', full.names = T)
ls.snp.names <- list.files('results/LP_DQ', 'SNP_step2.*regenie') %>% 
  gsub('.regenie|LP_SNP_step2_|_Oral|_NonOral|_All', '', .) 
ls.snp.names <- ls.snp.names %>% gsub('^A', 'A*', .) %>% gsub('^B', 'B*', .) %>% 
  gsub('_A', '_A*', .) %>% 
  gsub('_B', '_B*', .) %>% gsub('_DQB1', '_DQB1*', .) %>% gsub('_DQA1', '_DQA1*', .) %>% 
  gsub('_DRB1', '_DRB1*', .) %>% gsub('_', ' ', .)
ls.snp.names <- ls.snp.names %>% 
  gsub('^Class2', 'Class II:\nDQA1*0105 DQB1*0501 DQB1*0602 DRB1*0901 DQB1*0202 DQB1*0301', .)
ls.snp.names <- ls.snp.names %>% 
  gsub('^Class1', 'Class I:\nA*0201 A*0301 B*0801 B*1302 B*0702', .)
ls.snp.names <- ls.snp.names %>% 
  gsub('A* Class2', 'A & Class2', ., fixed = T)
ls.snp.names <- ls.snp.names %>% 
  gsub('B* Class2', 'B & Class2', ., fixed = T)


ls.snp.phenos <- rep(NA, length(ls.snp))
ls.snp.phenos[grepl('NonOral', ls.snp)] <- 'NonOral'
ls.snp.phenos[grepl('All', ls.snp)] <- 'All'
ls.snp.phenos[grepl('_Oral', ls.snp)] <- 'Oral'

ls.dat <- map(1:length(ls.snp), function(i) {
  tmp <- fread2(ls.snp[i])
  data.frame(tmp,
             Pheno = ls.snp.phenos[i],
             Adjusted = ls.snp.names[i])
}) %>% do.call(rbind, .) #%>% filter(Pheno != 'All')
ls.dat$Adjusted  <- factor(ls.dat$Adjusted, levels = ls.dat$Adjusted %>% unique %>% .[c(7,8,1,6,3:5,2)])

ls.dat$Position <- ls.dat$GENPOS/1e6
ls.dat$FDR_BY <- p.adjust(ls.dat$LOG10P %>% convert.p, method = 'BY')
  
# arrange phenos by beta
ls.dat$Pheno <- factor(ls.dat$Pheno, levels=c('Oral', 'All', 'NonOral'))

# positions of HLA genes
hla.pos <- fread2('data/HLA_gene_pos.tsv') %>% arrange(start) %>% .[-9, ]
hla.pos$start <- hla.pos$start/1e6
hla.pos$gene[8] <- 'HLA-DP'

# panel annotation
panel.anno <- data.frame(Position = 29.20, LOG10P = 124.5, Pheno = ls.dat$Pheno[1],
                         label = c('Class I', 'Class II', 'A & Class II', 'B & class II',
                                   'Class I & II without B*0702, DRB1*0901 and DQB1*0202',
                                   'Class I & II without B*0702 and DQB1*0202',
                                   'Class I & II without B*0702',
                                   'Class I and II'),
                         Adjusted = ls.dat %>% .$Adjusted %>% levels() %>% factor) 


# plot
p.hla.snp.manh <- ls.dat %>% filter(LOG10P > 2) %>% ggplot(aes(Position, LOG10P, color = Pheno)) +
  geom_point(alpha = 1, shape =1, size = .7) +
  xlab('Chr6 position (Mb)') +
  ylab(expression("-log"[10]*"(p)")) +
  scale_color_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev,
                     labels = c('non-OLP', 'All', 'OLP') %>% rev) + 
  scale_y_continuous(transform = 'log10', limits = c(2, 150), expand = c(0, 0.05)) +
  scale_x_continuous(breaks = c(seq(29, 34), hla.pos$start),
                     labels = c(seq(29, 34), paste0(c('\n', '\n', '\n\n', '\n', '\n\n', '\n\n\n', 
                                                      '\n\n\n\n', '\n\n\n\n\n'), 
                                                    hla.pos$gene))) +
  geom_vline(xintercept = hla.pos$start, color = 'grey50', linewidth = 0.2, linetype = 'dotted') +
  geom_hline(yintercept = filter(ls.dat, FDR_BY<0.05)$LOG10P %>% min(),
             linetype = 'solid', color = 'red') +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 2))) +
  facet_wrap(~ Adjusted, ncol = 2, scales = 'free_y') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, angle=0, hjust=0.5),
        axis.text.y=element_text(size=10),
        panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        axis.line=element_blank(),
        legend.position='bottom',
        legend.title=element_blank(),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), 'cm'),
        panel.border=element_rect(fill=NA, linewidth=0.3)) 

p.hla.snp.manh <- p.hla.snp.manh + geom_text(data = panel.anno, aes(label = label), color = 'black', hjust = 0) 

## Figure 4
jpeg('results/LP_DQ/Fig4_SNP_manhattan.jpg', width=11, height=10, units='in', res=1000)
p.hla.snp.manh %>% tag_facet(open = '', size = 4.3)
dev.off()
