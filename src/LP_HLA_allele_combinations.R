
## ===============================================
##
## Lichen planus HLA fine-mapping in R12 
## HLA allele combination analyses
##
## ===============================================

source('src/LP_DQ_functions.R')


## -----------------------------------------------
## DQA1-DQB1 diploid haplotypes 
## -----------------------------------------------

# read DQ pairs
dq <- fread('data/DQ.tsv') %>% unite(., Pair, A1, B1)
dq$Pair <- gsub('\\:', '', dq$Pair)

# allele dosages
h1 <- hla[, c(2, which(grepl('DQA1|DQB1', colnames(hla))))]
h1[, -1][ h1[, -1] < 0.5 ] <- 0
h1[, -1][ h1[, -1] > 1.5 ] <- 2
h1[, -1][ h1[, -1] >= 0.5 & h1[, -1] < 1.5  ] <- 1

# define DQA1-DQB1 allele pairs
dqb1.pairs <- expand.grid(colnames(h1)[grepl('DQA1', colnames(h1))], 
                          colnames(h1)[grepl('DQB1', colnames(h1))],
                          colnames(h1)[grepl('DQA1', colnames(h1))], 
                          colnames(h1)[grepl('DQB1', colnames(h1))])
dqb1.pairs <- unite(dqb1.pairs, PAIR1, Var1, Var2, remove=F)
dqb1.pairs <- unite(dqb1.pairs, PAIR2, Var3, Var4, remove=F)
dqb1.pairs <- unite(dqb1.pairs, PAIR, Var1, Var2, Var3, Var4, remove=F)
dqb1.pairs$PAIR <- dqb1.pairs$PAIR %>% gsub('\\:|', '', .) %>% gsub('\\*', '.', .)
dqb1.pairs$PAIR1 <- dqb1.pairs$PAIR1 %>% gsub('\\:|', '', .) %>% gsub('\\*', '.', .) %>% 
  gsub('DQA1.|DQB1.', '', .)
dqb1.pairs$PAIR2 <- dqb1.pairs$PAIR2 %>% gsub('\\:|', '', .) %>% gsub('\\*', '.', .) %>% 
  gsub('DQA1.|DQB1.', '', .)
dqb1.pairs <- dqb1.pairs %>% filter(PAIR1 %in% dq$Pair & PAIR2 %in% dq$Pair) %>% data.frame

# DQ pair dosages
dqb1.pair.geno <- map_dfc(1:nrow(dqb1.pairs), function(i) {
  # i <- 81
  alleles <- dqb1.pairs[i, c('Var1', 'Var2', 'Var3', 'Var4')] %>% unlist %>% as.character
  # alleles %>% unique %>% length
  out <- h1[, alleles]
  
  if((alleles %>% unique %>% length)==2) {
    out$TMP <- (out[, 1]==2 & out[, 2]==2)+0
  }
  if((alleles %>% unique %>% length)==3) {
    rep.allele <- colnames(out)[grepl('\\.1', colnames(out))] %>% 
      str_split_fixed(., '\\.', 2) %>% .[1]
    ind <- which(colnames(out)[1:4]==rep.allele)
    out <- out[, c(ind, c(1:4)[-ind])]
    ind <- which(grepl('\\.1', colnames(out)))
    out <- out[, -ind]
    out$TMP <- (out[, 1] == 2 & out[, 2] == 1 & out[, 3] == 1)+0
  }
  if((alleles %>% unique %>% length)==4) {
    out$TMP <- (out[, 1] == 1 & out[, 2] == 1 & out[, 3] == 1 & out[, 4] == 1)+0
  }
  
  colnames(out)[ncol(out)] <- dqb1.pairs[i, 'PAIR']
  return(out[, ncol(out)])
})

colnames(dqb1.pair.geno) <- dqb1.pairs$PAIR
asums <- dqb1.pair.geno %>% colSums(na.rm=T)
asums %>% hist(50)
which(asums<100)
dqb1.pair.geno <- dqb1.pair.geno[, asums>100]
colnames(dqb1.pair.geno) <- gsub('DQA1.', 'DQA1', colnames(dqb1.pair.geno)) %>% 
  gsub('DQB1.', 'DQB1', .)

fwrite(data.frame(IID = h1$IID, dqb1.pair.geno), 'data/genotypes/DQA1_DQB1_diplos.tsv', sep='\t')

# calculate n
dqb1.pair.geno %>% head
phenos %>% head
dqb1.pair.phenos <- left_join(phenos[, -1], data.frame(IID=h1[, 1], dqb1.pair.geno))
dqb1.pair.phenos <- data.frame(
  DQA1_DQB1_Pair = dqb1.pair.phenos %>% .[, -c(1:4)] %>% colnames,
  All = dqb1.pair.phenos %>% filter(All==1) %>% .[, -c(1:4)] %>% 
    apply(., 2, function(x) sum(x[x>0], na.rm=T)),
  NonOral = dqb1.pair.phenos %>% filter(NonOral==1) %>% .[, -c(1:4)] %>% 
    apply(., 2, function(x) sum(x[x>0], na.rm=T)),
  Oral = dqb1.pair.phenos %>% filter(Oral==1) %>% .[, -c(1:4)] %>% 
    apply(., 2, function(x) sum(x[x>0], na.rm=T))
)
fwrite(dqb1.pair.phenos, 'data/phenotypes/LP_DQA1_DQB1_diplos_phenos_n.tsv', sep='\t')

# convert to plink genotype format
createPlinkGeno(data.frame(IID=h1$IID, dqb1.pair.geno), 'LP_DQA1_DQB1_diplos')





## -----------------------------------------------
## regenie step 1: null model
## v3.0.1
## -----------------------------------------------

# Null model with basic covariates and all 3 phenotypes
system(paste0("regenie3 --step 1 --bed /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12_hm3 ",
              " --extract data/LP_regenie_null_qc_pass.snplist ", 
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars ",
              " --bt --bsize 1000 --lowmem --lowmem-prefix tmp/regenie_tmp_preds ",
              " --out results/LP_DQ/hm3_step1_Basic_All"))


## -----------------------------------------------
## regenie step 2: association 
## v3.0.1
## -----------------------------------------------


# DQA1-DQB1 diploid haplos
system(paste0("regenie3 --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_Basic_All_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQA1_DQB1diplos"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_Class2_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1 ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQA1_DQB1diplos_Class1"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_Class1_DQB10501hom_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10501hom ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_Class1_DQB10501hom"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_DQB10501hom_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10501hom ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_DQB10501hom"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_Class1_DQB10501het_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10501het ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_Class1_DQB10501het"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_DQB10501het_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10501het ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_DQB10501het"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_Class1_DQB10602hom_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10602hom ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_Class1_DQB10602hom"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_DQB10602hom_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10602hom ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_DQB10602hom"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_Class1_DQB10602het_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.Class1_DQB10602het ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_Class1_DQB10602het"))

system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bed data/genotypes/LP_DQA1_DQB1_diplos_haplos ",
              " --firth --approx --pThresh 0.01 ",
              " --pred results/LP_DQ/hm3_step1_DQB10602het_pred.list ",
              " --phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars.DQB10602het ",
              " --bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              " --out results/LP_DQ/LP_step2_DQB1diplos_DQB10602het"))



## ----------------------------------------------
## frequencies
## ----------------------------------------------

# LP controls
ind.controls <- phenos %>% filter(All == 0) %>% .$IID
ind.controls <- which(hla$IID %in% ind.controls)

sqrt(sum(hla[ind.controls, 'DQA1*01:01']>1.5 & 
           hla[ind.controls, 'DQB1*05:01']>1.5) / nrow(hla[ind.controls,]))

sqrt(sum(hla[ind.controls, 'DQA1*01:05']>1.5 & 
           hla[ind.controls, 'DQB1*05:01']>1.5) / nrow(hla[ind.controls,]))

sqrt(sum(hla[ind.controls, 'DRB1*15:01']>1.5 & 
           hla[ind.controls, 'DQA1*01:02']>1.5 & 
           hla[ind.controls, 'DQB1*06:02']>1.5) / nrow(hla[ind.controls, ]))


## ----------------------------------------------
## plot DQA1-DQB1 diploid haplo pairs
## ----------------------------------------------

# n
dqa1.dqb1.diplo.phenos.n <- fread('data/phenotypes/LP_DQA1_DQB1_diplos_phenos_n.tsv')
dqa1.dqb1.diplo.phenos.n$DQA1_DQB1_Pair <- dqa1.dqb1.diplo.phenos.n$DQA1_DQB1_Pair %>% 
  gsub('_B', ' / B', ., fixed=T) %>% 
  gsub('_A', ' | A', ., fixed=T) %>% 
  gsub('A1', 'DQA1*', ., fixed=T) %>% 
  gsub('B1', 'DQB1*', ., fixed=T) 
dqa1.dqb1.diplo.phenos.n$DQA1_DQB1_Pair <- gsub('/', '~', dqa1.dqb1.diplo.phenos.n$DQA1_DQB1_Pair,
                                                fixed = T) %>% gsub('|', '/', ., fixed = T)

# assoc data
res.dqa1.dqb1.diplos <- 
  map2(list.files('results/LP_DQ', 'LP_step2_DQA1_DQB1diplos_Class1_[A|O|N]+', full.names=T),
       c('All', 'NonOral', 'Oral'), function(d, p) process_result(d, p)) %>% 
  do.call(rbind, .)
res.dqa1.dqb1.diplos$ID <- res.dqa1.dqb1.diplos$ID %>% 
  gsub('_B', ' / B', ., fixed=T) %>% 
  gsub('_A', ' | A', ., fixed=T) %>% 
  gsub('A1', 'DQA1*', ., fixed=T) %>% 
  gsub('B1', 'DQB1*', ., fixed=T) 
res.dqa1.dqb1.diplos$ID <- gsub('/', '~', res.dqa1.dqb1.diplos$ID, fixed = T) %>% gsub('|', '/', ., fixed = T)

# arrange phenos
res.dqa1.dqb1.diplos$Pheno <- factor(res.dqa1.dqb1.diplos$Pheno, levels=c('Oral', 'All', 'NonOral'))

# select allele diplos in which at least one pheno is significant
res.dqa1.dqb1.diplos$Signif <- ifelse(res.dqa1.dqb1.diplos$FDR_BY<0.05, 'FDR <0.05', 'not significant')
res.dqa1.dqb1.diplos <- res.dqa1.dqb1.diplos %>% filter(ID %in% (res.dqa1.dqb1.diplos %>% group_by(ID) %>% 
                                                                   summarise(diploSignif=any(Signif=='FDR <0.05')) %>%
                                                                   .$ID %>% as.character))
# filter(diploSignif==T) %>% .$ID %>% as.character))
# include n
res.dqa1.dqb1.diplos <- left_join(res.dqa1.dqb1.diplos, dqa1.dqb1.diplo.phenos.n, by=c('ID'='DQA1_DQB1_Pair'))

# arrange by beta
res.dqa1.dqb1.diplos$ID <- factor(res.dqa1.dqb1.diplos$ID, levels=res.dqa1.dqb1.diplos %>% group_by(ID) %>% 
                                    summarise(BMean=mean(BETA)) %>% arrange(BMean) %>% .$ID)

# remove duplicates
res.dqa1.dqb1.diplos <- res.dqa1.dqb1.diplos %>% filter(!duplicated(UPPER))

# arrange ID label

# create data frames for plotting
res.dqa1.dqb1.diplos.dqb1.0501 <- res.dqa1.dqb1.diplos %>% filter(grepl('DQB1.0501', ID)) %>% 
  filter(!grepl('DQB1.0602', ID), All > 30)
diplo.parts <- res.dqa1.dqb1.diplos.dqb1.0501$ID %>% str_split_fixed(., ' \\| ', 2)
res.dqa1.dqb1.diplos.dqb1.0501$ID[grepl('DQB1.0501', diplo.parts[, 2])] <- 
  diplo.parts[grepl('DQB1.0501', diplo.parts[, 2]), 2:1] %>% data.frame %>% unite(., U, sep = ' | ') %>% .$U
res.dqa1.dqb1.diplos.dqb1.0501.n <- res.dqa1.dqb1.diplos.dqb1.0501 %>% 
  pivot_longer(cols=c('All', 'NonOral', 'Oral'))
res.dqa1.dqb1.diplos.dqb1.0501.n$name <- factor(res.dqa1.dqb1.diplos.dqb1.0501.n$name, 
                                                levels=c('Oral', 'All', 'NonOral'))

res.dqa1.dqb1.diplos.dqb1.0602 <- res.dqa1.dqb1.diplos %>% filter(grepl('DQB1.0602', ID), All > 30)
diplo.parts <- res.dqa1.dqb1.diplos.dqb1.0602$ID %>% str_split_fixed(., ' \\| ', 2)
res.dqa1.dqb1.diplos.dqb1.0602$ID[grepl('DQB1.0602', diplo.parts[, 2])] <- 
  diplo.parts[grepl('DQB1.0602', diplo.parts[, 2]), 2:1] %>% data.frame %>% unite(., U, sep = ' | ') %>% .$U
res.dqa1.dqb1.diplos.dqb1.0602.n <- res.dqa1.dqb1.diplos.dqb1.0602 %>% 
  pivot_longer(cols=c('All', 'NonOral', 'Oral'))
res.dqa1.dqb1.diplos.dqb1.0602.n$name <- factor(res.dqa1.dqb1.diplos.dqb1.0602.n$name, 
                                                levels=c('Oral', 'All', 'NonOral'))

res.dqa1.dqb1.diplos.other   <- res.dqa1.dqb1.diplos %>% filter(!grepl('DQB1.0501|DQB1.0602', ID))
res.dqa1.dqb1.diplos.other.n <- res.dqa1.dqb1.diplos.other %>% 
  pivot_longer(cols=c('All', 'NonOral', 'Oral'))
res.dqa1.dqb1.diplos.other.n$name <- factor(res.dqa1.dqb1.diplos.other.n$name, 
                                            levels=c('Oral', 'All', 'NonOral'))

# legend
p.legend <- res.dqa1.dqb1.diplos.dqb1.0501 %>% 
  ggplot(aes(OR, ID, xmin=OR_LOWER, xmax=OR_UPPER, color=Pheno, shape=Signif)) +
  geom_pointrange(position=position_dodge(width=0.6), size=0.6, linewidth=0.7) +
  xlab('beta') + ylab('') +
  scale_color_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev,
                     labels = c('non-OLP', 'All', 'OLP') %>% rev) + 
  scale_shape_manual(name=NULL, values=c(15, 0)) +
  geom_vline(xintercept=0, linewidth=0.3, color='black', linetype='dashed') +
  guides(color = guide_legend(override.aes=list(shape = 15))) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10, margin=margin(1,3,1,1)),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3),
        legend.position='bottom',
        legend.title=element_blank(),
        plot.margin=unit(c(.5,0.0,.5,.5), 'cm'),
        panel.border=element_rect(fill=NA, linewidth=0.2))
p.legend <- p.legend %>% get_legend() %>% as_ggplot()

# plot 0501
p.res.dqa1.dqb1.diplos.dqb1.0501 <- res.dqa1.dqb1.diplos.dqb1.0501 %>% 
  ggplot(aes(OR, ID, xmin=OR_LOWER, xmax=OR_UPPER, color=Pheno, shape=Signif)) +
  geom_pointrange(position=position_dodge(width=0.6), size=0.6, linewidth=0.7) +
  xlab('OR') + ylab('') +
  scale_color_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev,
                     labels = c('non-OLP', 'All', 'OLP') %>% rev) + 
  scale_shape_manual(name=NULL, values=c(15, 0)) +
  scale_x_continuous(transform = 'log10') +
  geom_vline(xintercept=1, linewidth=0.3, color='black', linetype='dashed') +
  guides(color = guide_legend(override.aes=list(shape = 15))) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=8.4, margin=margin(1,3,1,1)),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3),
        legend.position='none',
        legend.title=element_blank(),
        plot.margin=unit(c(.5,0.0,.5,.5), 'cm'),
        panel.border=element_rect(fill=NA, linewidth=0.2))

sapply(seq(1, unique(res.dqa1.dqb1.diplos.dqb1.0501$ID) %>% length, by=2), function(i) {
  p.res.dqa1.dqb1.diplos.dqb1.0501 <<- p.res.dqa1.dqb1.diplos.dqb1.0501 +
    annotate(geom='rect', ymin=i-0.5, ymax=i+0.5, xmin= 0, xmax=Inf,
             fill='grey25', alpha=0.05, color=NA)
})

p.res.dqa1.dqb1.diplos.dqb1.0501.n <- res.dqa1.dqb1.diplos.dqb1.0501.n %>% 
  ggplot(aes(value, ID, fill=name)) +
  geom_bar(stat='identity', position=position_dodge(width=0.6), color='NA', width=0.5) +
  scale_fill_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev,
                    labels = c('non-OLP', 'All', 'OLP') %>% rev) + 
  xlab('# cases') + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid=element_blank(),
        legend.position='none',
        plot.margin=unit(c(.5,.5,.5,0.0), 'cm'))

sapply(seq(1, unique(res.dqa1.dqb1.diplos.dqb1.0501.n$ID) %>% length, by=2), function(i) {
  p.res.dqa1.dqb1.diplos.dqb1.0501.n <<- p.res.dqa1.dqb1.diplos.dqb1.0501.n +
    annotate(geom='rect', ymin=i-0.5, ymax=i+0.5, xmin= -Inf, xmax=Inf,
             fill='grey25', alpha=0.05, color=NA)
})

p.res.dqa1.dqb1.diplos.dqb1.0501 + p.res.dqa1.dqb1.diplos.dqb1.0501.n + 
  plot_layout(widths=c(2.8, 1)) 


# plot 0602
p.res.dqa1.dqb1.diplos.dqb1.0602 <- res.dqa1.dqb1.diplos.dqb1.0602 %>% 
  ggplot(aes(OR, ID, xmin=OR_LOWER, xmax=OR_UPPER, color=Pheno, shape=Signif)) +
  geom_pointrange(position=position_dodge(width=0.6), size=0.6, linewidth=0.7) +
  xlab('OR') + ylab('') +
  scale_color_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev,
                     labels = c('non-OLP', 'All', 'OLP') %>% rev) + 
  scale_shape_manual(name=NULL, values=c(15, 0)) +
  scale_x_continuous(trans = 'log10') +
  geom_vline(xintercept=1, linewidth=0.3, color='black', linetype='dashed') +
  guides(color = guide_legend(override.aes=list(shape = 15))) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=8.4, margin=margin(1,3,1,1)),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3),
        legend.position='none',
        legend.title=element_blank(),
        plot.margin=unit(c(.5,0.0,.5,.5), 'cm'),
        panel.border=element_rect(fill=NA, linewidth=0.2))

sapply(seq(1, unique(res.dqa1.dqb1.diplos.dqb1.0602$ID) %>% length, by=2), function(i) {
  p.res.dqa1.dqb1.diplos.dqb1.0602 <<- p.res.dqa1.dqb1.diplos.dqb1.0602 +
    annotate(geom='rect', ymin=i-0.5, ymax=i+0.5, xmin= 0, xmax=Inf,
             fill='grey25', alpha=0.05, color=NA)
})

p.res.dqa1.dqb1.diplos.dqb1.0602.n <- res.dqa1.dqb1.diplos.dqb1.0602.n %>% 
  ggplot(aes(value, ID, fill=name)) +
  geom_bar(stat='identity', position=position_dodge(width=0.6), color='NA', width=0.5) +
  scale_fill_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev,
                    labels = c('non-OLP', 'All', 'OLP') %>% rev) + 
  xlab('# cases') + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid=element_blank(),
        legend.position='none',
        plot.margin=unit(c(.5,.5,.5,0.0), 'cm'))

sapply(seq(1, unique(res.dqa1.dqb1.diplos.dqb1.0602.n$ID) %>% length, by=2), function(i) {
  p.res.dqa1.dqb1.diplos.dqb1.0602.n <<- p.res.dqa1.dqb1.diplos.dqb1.0602.n +
    annotate(geom='rect', ymin=i-0.5, ymax=i+0.5, xmin= 0, xmax=Inf,
             fill='grey25', alpha=0.05, color=NA)
})

p.res.dqa1.dqb1.diplos.dqb1.0602 + p.res.dqa1.dqb1.diplos.dqb1.0602.n + 
  plot_layout(widths=c(2.8, 1)) 


# combined plot
p.res.dqa1.dqb1.diplos.2 <- ggpubr::ggarrange(
  ggpubr::ggarrange(p.res.dqa1.dqb1.diplos.dqb1.0501 + p.res.dqa1.dqb1.diplos.dqb1.0501.n +
                      plot_layout(widths=c(2.5, 1)),
                    (p.res.dqa1.dqb1.diplos.dqb1.0602 + p.res.dqa1.dqb1.diplos.dqb1.0602.n +
                       plot_layout(widths=c(2.5, 1))),
                    nrow = 2, labels = c('a)', 'b)'), font.label = list(size = 14), heights = c(45, 36)),
  p.legend, nrow = 2, heights = c(1, 0.05))


## Figure 2

jpeg('results/LP_DQ/Fig2_DQA1_DQB1_diploid_haplos.jpg', width=8, height=11, units='in', res=1000)
p.res.dqa1.dqb1.diplos.2
dev.off()

