
## ===============================================
##
## Lichen planus HLA fine-mapping in R12 
## HLA allele homozygote analyses
##
## ===============================================

source('src/LP_DQ_functions.R')


## -----------------------------------------------
## List of samples that are DQA1-DQB1 homozygotes
## -----------------------------------------------

# list all possible DQA1-DQB1 pairs
dqa1.dqb1.pairs <- expand.grid(
  dplyr::select(hla, contains('DQA1')) %>% colnames,
  dplyr::select(hla, contains('DQB1')) %>% colnames
)
colnames(dqa1.dqb1.pairs) <- c('DQA1', 'DQB1')
dqa1.dqb1.pairs <- unite(dqa1.dqb1.pairs, 'Pair', sep='_', remove=F)

# list of samples for each DQ pair homozygote
dqa1.dqb1.pairs.samples <- map(1:nrow(dqa1.dqb1.pairs), function(i) {
  hla[
    hla[, as.character(dqa1.dqb1.pairs[i, 'DQA1'])] > 1.5 &
      hla[, as.character(dqa1.dqb1.pairs[i, 'DQB1'])] > 1.5,  'IID']
})
names(dqa1.dqb1.pairs.samples) <- dqa1.dqb1.pairs$Pair

# number of LP cases for each homozygote pair
dqa1.dqb1.pairs.cases <- map(1:length(dqa1.dqb1.pairs.samples), function(i) {
  intersect(dqa1.dqb1.pairs.samples[[i]], filter(phenos, All==1)$IID) %>% length
}) %>% unlist
dqa1.dqb1.pairs <- data.frame(dqa1.dqb1.pairs, Cases=dqa1.dqb1.pairs.cases)
# keep pairs with at least some samples
dqa1.dqb1.pairs <- dqa1.dqb1.pairs %>% filter(Cases>5)

# updated list of samples for each DQ homozygote
dqa1.dqb1.pairs.samples <- map(1:nrow(dqa1.dqb1.pairs), function(i) {
  hla[ 
    hla[, as.character(dqa1.dqb1.pairs[i, 'DQA1'])] > 1.5 & 
      hla[, as.character(dqa1.dqb1.pairs[i, 'DQB1'])] > 1.5,  'IID']
})
names(dqa1.dqb1.pairs.samples) <- dqa1.dqb1.pairs$Pair

# number of samples in each hom group
dqa1.dqb1.pairs.samples %>% map(length) %>% data.frame %>% t %>% 
  data.frame %>% rownames_to_column() %>% rename('haplo'='rowname', 'n'='.') %>% 
  fwrite('data/phenotypes/DQhomozygote_samples_n.tsv', sep='\t')
fread('data/phenotypes/DQhomozygote_samples_n.tsv')

# number of cases 
dqa1.dqb1.pairs.samples.phenos <-  map2_dfr(dqa1.dqb1.pairs.samples, 
          dqa1.dqb1.pairs.samples %>% names, function(x, y) {
  data.frame(
    Haplo=y,
    All=intersect(filter(phenos, All==1)$IID, x) %>% length,
    Oral=intersect(filter(phenos, Oral==1)$IID, x) %>% length,
    NonOral=intersect(filter(phenos, NonOral==1)$IID, x) %>% length)
})
dqa1.dqb1.pairs.samples.phenos$Haplo <- gsub('\\:', '', dqa1.dqb1.pairs.samples.phenos$Haplo) %>%  gsub('_', ' | ', .)
fwrite(dqa1.dqb1.pairs.samples.phenos, 'data/phenotypes/dqa1.dqb1.pairs.samples.phenos', sep='\t')
fread('data/phenotypes/dqa1.dqb1.pairs.samples.phenos') %>% arrange(1/All) %>% kable(format='simple')

# write the DQ homozygote set of samples
data.frame(FID=dqa1.dqb1.pairs.samples %>% unlist %>% unique,
           IID=dqa1.dqb1.pairs.samples %>% unlist %>% unique) %>% 
  fwrite(., 'data/phenotypes/DQhomozygote_samples.list', sep='\t')
fread( 'data/phenotypes/DQhomozygote_samples.list')

# removing effect hz's
tmp <- fread('data/phenotypes/DQhomozygote_samples.list') %>% 
  filter(!(IID %in% c(dqa1.dqb1.pairs.samples[['DQA1*01:01_DQB1*05:01']], 
                      dqa1.dqb1.pairs.samples[['DQA1*01:02_DQB1*06:02']])))
fwrite(tmp, 'data/phenotypes/DQhomozygotes_samples_noeff.list', sep='\t')



## -----------------------------------------------
## regenie step 2: association
## v3.0.1
## -----------------------------------------------

# basic = HLA association without any additional allele covariates
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_Basic_All_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno --covarFile data/phenotypes/LP.covars ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_basic"))

# HLA association adjusted for DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_DQB10501"))

# HLA association adjusted for DQA1*0101, DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10101_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10101_DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_DQA10101_DQB10501"))

# HLA association adjusted for DQA1*0105, DQB1*0501 dosages 
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_DQA10105_DQB10501"))

# association adjusted for DQA1*0105, DQB1*0602 and DQB1*0501
system(paste0("/home/ivm/Documents/tmp/regenie_v3.6.gz_x86_64_Linux ",
              " --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_DQA10105_DQB10501_DQB10602"))

# association adjusted for DQA1*0105, DQA1*0101, DQB1*0501
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQB10501_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_DQA10105_DQA10101_DQB10501"))

# association adjusted for DQA1*0105, DQA1*0101, DQB1*0501, DQB1*0602
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_DQA10105_DQA10101_DQB10501_DQB10602"))

# association adjusted for DQA1*0105, DQA1*0101, DQA1*0102, DQB1*0501, DQB1*0602
system(paste0("regenie3 --step 2 --bgen /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen ",
              "--sample /finngen/library-red/finngen_R12/hla_1.0/bgen/R12_HLA.bgen.sample ",
              "--keep data/phenotypes/DQhomozygote_samples.list ",
              "--ref-first --firth --approx --pThresh 0.0001 ",
              "--pred results/LP_DQ/hm3_step1_DQA10105_DQA10101_DQA10102_DQB10501_DQB10602_pred.list ",
              "--phenoFile data/phenotypes/LP.pheno ",
              "--covarFile data/phenotypes/LP.covars.DQA10105_DQA10101_DQA10102_DQB10501_DQB10602 ",
              "--bt --bsize 400 --lowmem --lowmem-prefix tmp/regenie_tmp_step2 ",
              "--out results/LP_DQ/LP_HLA_DQhom_step2_DQA10105_DQA10101_DQA10102_DQB10501_DQB10602"))



## -----------------------------------------------
## results
## -----------------------------------------------

# list regenie assoc data
res <- list(list.files('results/LP_DQ', 'LP_HLA_DQhom_step2_basic_[A|O|N]+', full.names=T),
            #list.files('results/LP_DQ', 'LP_HLA_DQhom_step2_DQB10501_[A|O|N]+', full.names=T),
            #list.files('results/LP_DQ', 'LP_HLA_DQhom_step2_DQA10101_DQB10501_[A|O|N]+', full.names=T),
             list.files('results/LP_DQ', 'LP_HLA_DQhom_step2_DQA10105_DQB10501_[A|O|N]+', full.names=T),
            list.files('results/LP_DQ', 'LP_HLA_DQhom_step2_DQA10105_DQB10501_DQB10602_[A|O|N]+', full.names=T))

# adjustment alleles in the same order
adjustment.chain <- c('none', 'DQA1*0105 & DQB1*0501', 'DQA1*0105 & DQB1*0501 &\nDQB1*0602')
                   
target.alleles <- c('A*03:01', 'B*13:02', 'B*35:01', 'DQA1*01:02', 'B*08:01', 'DQB1*03:01',
                    'DQB1*02:01', 'DQB1*05:01', 'DQB1*02:02', 'DQB1*06:02', 'DRB1*09:01', 'B*07:02')

# read allele results data
res.alleles <- map_dfr(1:length(res), function(i) {
  # i <- 4
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
#res.alleles %>% group_by(Pheno) %>% summarise(MB=mean(BETA)) %>% arrange(MB) %>% .$Pheno %>% rev)

panel.anno <- data.frame(OR = 0.55, ID = 13.25, Pheno = res.alleles$Pheno[1], OR_LOWER = 2, OR_UPPER = 2, 
                         Signif = res.alleles$Signif[1],
                         label = rep(paste0(letters[1:3], ')')),
                         Adjusted = res.alleles$Adjusted %>% as.character %>% unique %>% factor) 


## -----------------------------------------------
## plots
## -----------------------------------------------

p.res.alleles.DQhom <- res.alleles %>% 
  ggplot(aes(OR, ID, xmin=OR_LOWER, xmax=OR_UPPER, color=Pheno, shape=Signif)) +
  geom_pointrange(position=position_dodge(width=0.6), size=0.6, linewidth=0.7) +
  xlab('OR') + ylab('') +
  scale_color_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev) + # #FFD580 ff9248
  scale_shape_manual(name=NULL, values=c(15, 0)) +
  scale_x_continuous(transform = 'log10') +
  coord_cartesian(xlim = c(0.55, 2.5), ylim = c(1, 12), clip = 'off') +
  geom_text(data = panel.anno, aes(OR, ID, label = label, xmin = OR_LOWER, xmax = OR_UPPER, shape = Signif), 
            size = 4.3, color = 'black', fontface = 'bold') +
  geom_vline(xintercept = 1, linewidth = 0.3, color = 'black', linetype = 'dashed') +
  guides(color = guide_legend(override.aes=list(shape = 15))) +
  facet_wrap(vars(Adjusted), scales='fixed', ncol=3) +
  theme_minimal() +
  #scale_x_continuous(limits = c(-0.45, 0.69), oob = scales::squish()) +
  theme(axis.text.x=element_text(size=10, angle=0, hjust=0),
        axis.text.y=element_text(size=10),
        panel.grid=element_blank(),
        axis.line=element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3, colour = 'black'),
        legend.position='bottom',
        legend.title=element_blank(),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), 'cm'),
        panel.border=element_rect(fill=NA, linewidth=0.2))
sapply(seq(1, unique(res.alleles$ID) %>% length, by=2), function(i) {
  p.res.alleles.DQhom <<- p.res.alleles.DQhom +
    annotate(geom='rect', ymin=i-0.5, ymax=i+0.5, xmin= 0, xmax=Inf, 
             fill='grey25', alpha=0.05, color=NA) 
}) 

## number of cases in hom haplotypes
dqa1.dqb1.pairs.samples.phenos <- fread('data/phenotypes/dqa1.dqb1.pairs.samples.phenos') %>% 
  pivot_longer(2:4)
dqa1.dqb1.pairs.samples.phenos$Haplo <- factor(dqa1.dqb1.pairs.samples.phenos$Haplo,
                                               levels=dqa1.dqb1.pairs.samples.phenos %>% group_by(Haplo) %>% 
                                                 summarise(MM=mean(value)) %>% arrange(MM) %>% .$Haplo)
dqa1.dqb1.pairs.samples.phenos$name <- factor(dqa1.dqb1.pairs.samples.phenos$name, levels=c('Oral', 'All', 'NonOral'))

p.dqa1.dqb1.pairs.samples.phenos <- dqa1.dqb1.pairs.samples.phenos %>% ggplot(aes(value, Haplo, fill=name)) +
  geom_bar(stat='identity', position=position_dodge(width=0.6), color='NA', width=0.5) +
  scale_fill_manual(values=c("#440154FF", "#21908CFF", "#FFBF00") %>% rev) +
  scale_x_continuous(transform = 'log10') +
  coord_cartesian(clip = 'off', ylim = c(0.5, 10.5), xlim = c(1, 1200), expand = F) +
  xlab('# cases') + 
  geom_text(x = 0.065, y = 11.32, label = 'd)', size = 4.3, color = 'black', fontface = 'bold') +
  ylab('') +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10, margin=margin(1,3,1,1)),
        axis.title.y=element_text(),
        panel.grid=element_blank(),
        axis.ticks.x = element_line(linewidth = 0.3, colour = 'black'),
        legend.position='right',
        #legend.position.inside=c(0.7, 0.3),
        #legend.background=element_rect(fill=NA, linewidth=0.2, color='grey'),
        legend.title=element_blank(),
        plot.margin=unit(c(.5,.5,.5, 0.5), 'cm'),
        panel.border=element_rect(fill=NA, linewidth=0.2))
sapply(seq(1, unique(dqa1.dqb1.pairs.samples.phenos$Haplo) %>% length, by=2), function(i) {
  p.dqa1.dqb1.pairs.samples.phenos <<- p.dqa1.dqb1.pairs.samples.phenos +
    annotate(geom='rect', ymin=i-0.5, ymax=i+0.5, xmin= 0, xmax=Inf, 
             fill='grey25', alpha=0.05, color=NA) 
}) 


## Figure 1
jpeg('results/LP_DQ/Fig1_Alleles_DQhom.jpg', width=11, height=8.6, units='in', res=1000)
p.res.alleles.DQhom / (p.dqa1.dqb1.pairs.samples.phenos + plot_spacer() + plot_layout(widths = c(1, 0.36))) + 
  plot_layout(heights = c(1, 0.65))
dev.off()







