# this will be the working code to process the 
# RNA-seq data from Xiao'ning who is from Jingzhou lab
# @author Yisong
# @since  2017-12-27
# @update 2018-01-19
#---

#---
# data QC was performed using shell script:
#  :
# the original data was transfered from mobile disk
# and figureprints were in consistant with the md5.txt
# from the sequencing company
#---

#---
# 1. study objectives
# 2. hypotheses
# 3. power and sample size calculation
# 4. data collection methods
# 5. validation methods
# 6. data summaries
# 7. detail of the statistical analysis
# 8. sources of missing data
# 9. exceptions
#---


#---
# 
# 1a. Ask a clinically relevant, testable question
# 2b. Choose an experimental design and statistical framework
# 3c. Set up a null hypothesis, that is, a testable claim 
#     that becomes the target of statistical analysis
# 4d. Fix a rejection region, that is, the degree of 
#     evidence against the null hypothesis at which it may
#     be rejected
# 5e. Conduct the experiment: collect
#     data, compute the test statistics
# 6f. Report results: all and only those
#     genes that fall within the rejection
#     region are declared differentially
#     expressed
#---


#---
# https://stackoverflow.com/questions/18306362/run-r-script-from-command-line
#---
pkgs <- c( 'tidyverse','Rsubread','org.Mm.eg.db','edgeR',
           'limma', 'DESeq2', 'genefilter','grid',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'cowplot', 'readxl', 'magrittr',
           'BSgenome.Mmusculus.UCSC.mm10',
           'BSgenome.Mmusculus.UCSC.mm10.Rbowtie')

load.lib <- lapply(pkgs, require, character.only = TRUE)

#---
# load processed data
#---

working.env <- 'window'
linux.path  <- file.path('/wa/zhenyisong/results/chenlab/xiaoning/data')
window.path <- file.path('D:\\yisong.data')


#---
# [zhenyisong@mgt data]$ md5sum xiaon.Rdata
# 954e0ee9a5969a873f2f0fd9791d0375  xiaon.Rdata
#---

switch( working.env, linux  = { setwd(linux.path);
                                load('xiaon.Rdata') },
                     window = { setwd(window.path);
                                load('xiaon.Rdata') } )


# statistical design
# 2 x 2 design
#---

design.table <- data.frame( 'LU'  = c('AL_LU', 'A_LU'),
                            'LU-' = c('AL', 'A'), check.names = FALSE)
rownames(design.table) <- c('LN','LN-')

design.tableGrob <- tableGrob( design.table)

grid.newpage()
grid.draw(design.tableGrob)


#---
# now process the raw data
#---


rsubread.index.lib   <- file.path('/home/zhenyisong/data/reference/index')

xiaon.file.path      <- file.path(
                          '/wa/zhenyisong/results/chenlab/xiaoning/data')
xiaon.output.dir     <- file.path(
                          '/home/zhenyisong/data/results/chenlab/xiaoning/qc.step')
xiaon.clean.data     <- list.files( path         = xiaon.file.path, 
                                    pattern      = '.fq.gz$', 
                                    all.files    = FALSE, 
                                    full.names   = TRUE, 
                                    recursive    = TRUE, 
                                    ignore.case  = FALSE, 
                                    include.dirs = TRUE)

read.1.files             <- grep('R1',xiaon.clean.data) %>% 
                            xiaon.clean.data[.]
read.2.files             <- grep('R2',xiaon.clean.data) %>% 
                            xiaon.clean.data[.]

xiaon.output.filenames   <- basename(read.1.files) %>% 
                            sub(pattern = '_R1.fq.gz', replacement = '') %>%
                            paste0(xiaon.output.dir,'/', . ,'.bam')
xioan.sample.names       <- basename(read.1.files) %>% 
                            sub(pattern = '_R1.fq.gz', replacement = '')

sampleFile               <- tempfile(pattern = 'zhen3temp', tmpdir = tempdir())
sample.file              <- data.frame( FileName1  = read.1.files,
                                        FileName2  = read.2.files, 
                                       SampleName = xioan.sample.names)

write_tsv(sample.file, path = sampleFile)
genome          <- 'BSgenome.Mmusculus.UCSC.mm10'
cluster         <- makeCluster(3)

setwd(xiaon.output.dir)

xioan.qPorject  <- qAlign( sampleFile,
                           genome,
                           auxiliaryFile = NULL,
                           aligner = 'Rbowtie',
                           maxHits = 1,
                           paired  = 'fr',
                           splicedAlignment = FALSE,
                           snpFile = NULL,
                           bisulfite = 'no',
                           alignmentParameter = NULL,
                           projectName = 'xioan',
                           alignmentsDir = xiaon.output.dir,
                           lib.loc  = NULL,
                           cacheDir = NULL,
                           clObj = cluster,
                           checkOnly = F)

qQCReport( xioan.qPorject, 
           pdfFilename    = 'xioan.QC.pdf', 
           useSampleNames = TRUE, 
           clObj          = cluster)
stopCluster(cluster)

#---
# module 2
# start rsubread to count the gene expression
#---

setwd(rsubread.index.lib)

base.string      <-  'mm10'

align.result     <- align( index          = base.string, 
                           readfile1      = read.1.files, 
                           readfile2      = read.2.files, 
                           input_format   = 'gzFASTQ', 
                           type           = 'rna',
                           output_file    = xiaon.output.filenames, 
                           output_format  = 'BAM',
                           PE_orientation = 'fr', 
                           nthreads       = 4, 
                           indels         = 1,
                           maxMismatches  = 3,
                           phredOffset    = 33,
                           unique         = T )
                           
xiaon.genes     <- featureCounts( xiaon.output.filenames, 
                                  useMetaFeatures        = TRUE,
                                  countMultiMappingReads = FALSE,
                                  strandSpecific         = 0, 
                                  isPairedEnd            = TRUE,
                                  requireBothEndsMapped  = TRUE,
                                  autosort               = TRUE,
                                  nthreads               = 4,
                                  annot.inbuilt          = 'mm10', 
                                  allowMultiOverlap      = TRUE)

#---
# this is the temprotmanet save the
# intermidiate results of the processing
# results
# I transferred the data to window disk in my laptop
# X1, yisong.data directory.
#---

'
setwd(xiaon.file.path)
save.image(file = 'xiaon.Rdata')
quit('no')
'

#---
#
# PCA post-QC check
# if the sample clustering 
# reach to the expectation
# or should we continue to proceed the data
# and perform differential expression analysis
#
#---


#colnames(xiaon.genes$counts)

sample.names   <- c( paste('AL_HU', '_',1:4, sep = ''),
                     paste('AL_ICH','_',1:5, sep = ''),
                     paste('AL_LU','_',1:4, sep = ''),
                     paste('A','_',1:4, sep = ''),
                     paste('A_LU','_',1:4, sep = ''),
                     paste('Veh','_',1:5, sep = '') )

mapping.QC             <- as.matrix(xiaon.genes$stat[,-1])
rownames(mapping.QC)   <- xiaon.genes$stat[,1]
colnames(mapping.QC)   <- sample.names
mapping.ratio.QC       <- {mapping.QC / colSums(mapping.QC)} %>%
                           format(scientific = T, digits = 3)
defined.theme          <- ttheme_default( base_size = 6)
mapping.QC.table       <- mapping.QC %>% 
                          tableGrob( theme = defined.theme )
mapping.ratio.QC.table <- tableGrob( mapping.ratio.QC, 
                                     theme = defined.theme)

grid.newpage()
grid.draw(mapping.QC.table)
grid.newpage()
grid.draw(mapping.ratio.QC.table)

samples.xiaon  <- data.frame( treat = c( rep('AL_HU',4),
                                         rep('AL_ICH',5),
                                         rep('AL_LU',4),
                                         rep('A', 4),
                                         rep('A_LU', 4),
                                         rep('Veh', 5)) )

xioan.vst     <- xiaon.genes %$% counts %>%
                 DESeqDataSetFromMatrix( colData = samples.xiaon, 
                                         design  = ~ treat) %>%
                 DESeq() %>%
                 varianceStabilizingTransformation() %>% 
                 assay()
colnames(xioan.vst) <- sample.names
sds <- rowSds(xioan.vst)
sh  <- shorth(sds)
(sh)
# [1] 0.2169877
xiaon.pca          <- xioan.vst %>% subset(sds > sh) %>%
                      t() %>% prcomp
xioan.pve          <- xiaon.pca$sdev^2/sum(xiaon.pca$sdev^2)
(xioan.pve)
xiaon.groups       <- sample.names

#---
# quick PCA plot for chen's work
#---
x1.range <- xiaon.pca$x[,1] %>% range
x2.range <- xiaon.pca$x[,2] %>% range
(x1.range) 
(x2.range)
plot( xiaon.pca$x[,c(1,2)], 
      xlim = c(-90,100), ylim = c(-70, 60), 
      col = 'orange', pch = 16 )
text( xiaon.pca$x[,1], xiaon.pca$x[,2], 
      labels = xiaon.groups , cex = 0.75, 
      pos    = 3, adj = c(0,1))

#---
# complicate ggplot2 map for PCA analysis
# see my manual in Hackseq project
#---
pca.no             <- dim(xiaon.pca$x)[2]
graph.text.size    <- 7
point.size         <- 1.5
scree.plot.ggplot <- {xiaon.pca$sdev^2} %>%
                     {./sum(.)} %>% 
                     {data.frame(variance = ., pca = c(1:pca.no)) } %>%
                     ggplot(data = .) +
                     xlab('PCA Scree Plot') +
                     ylab('Variance Proportion') +
                     scale_x_continuous( breaks = c(1:pca.no), 
                                         labels = as.character(c(1:pca.no), 
                                         limits = as.character(c(1:pca.no)))) +
                     geom_point(aes(x = pca, y = variance), size = point.size) +
                     geom_line(aes(x = pca, y = variance), size = 0.8) +
                     scale_linetype_discrete() +
                     theme(legend.position = 'none') +
                     theme_classic() +
                     theme( aspect.ratio = 1, text = element_text(size = graph.text.size))
(scree.plot.ggplot)

cumsum.plot.ggplot <- {xiaon.pca$sdev^2} %>%
                      {./sum(.)} %>% 
                      {data.frame(variance = cumsum(.), pca = c(1:pca.no))}%>%
                      ggplot(data = .) +
                      xlab('Principle Components') +
                      ylab('Culmulative Variance') +
                      scale_x_continuous( breaks = c(1:pca.no), labels = as.character(c(1:pca.no), 
                                          limits = as.character(c(1:pca.no)))) +
                      geom_point(aes(x = pca, y = variance), size = point.size) +
                      geom_line(aes(x = pca, y = variance),size = 0.8) +
                      scale_linetype_discrete() +
                      theme(legend.position = 'none') +
                      theme_classic() +
                      theme( aspect.ratio = 1, text = element_text(size = graph.text.size))
(cumsum.plot.ggplot)

xiaon.PCA.ggplot <- as.data.frame(xiaon.pca$x) %>% 
                    mutate(color.choice = { rownames(.) %>% 
                                            sub('_\\d+$', '', perl = F, .) %>%
                                            as.factor}) %>%
                    ggplot(data = . ) + 
                    geom_point(aes(x = PC1, y = PC2, color = color.choice), size = point.size) + 
                    scale_colour_manual( name   = 'group in experiment design',
                                         values = c( 'darkblue','brown3','yellow',
                                                     'deeppink','lemonchiffon3','green'),
                                         labels = c( 'A','A_LU','AL_HU',
                                                     'AL_ICH','AL_LU', 'Veh')) +
                    theme_classic() +
                    theme( aspect.ratio = 1, text = element_text(size = graph.text.size))
(xiaon.PCA.ggplot)

source.pca.data <-  as.data.frame(xiaon.pca$x) %>% 
                    mutate( xiaon.groups = { rownames(.) %>% 
                                             sub('_\\d+$', '', perl = F, .) %>%
                                             as.factor } )  %>%
                    as.tibble() %>%
                    dplyr::select(PC1,PC2,xiaon.groups)

source.pca.table <- grid.newpage() %>% 
                   {tableGrob( source.pca.data, rows = NULL)} %>%
                   grid.draw()


PCA.QC.final.ggplot <- ggdraw() +
                       draw_plot(xiaon.PCA.ggplot, x = 0, y = 0.2, width = 0.4, height = 0.4) +
                       draw_plot(scree.plot.ggplot, x = 0.4, y = 0.4, width = 0.2, height = 0.2) +
                       draw_plot(cumsum.plot.ggplot, x = 0.4, y = 0.2, width = 0.2, height = 0.2) +
                       draw_plot_label( label = LETTERS[1:3], 
                                        x = c(0, 0.4, 0.4), y = c(0.6, 0.6, 0.4),
                                        size = 9)
(PCA.QC.final.ggplot)
#---
# tissue specific genes
# checking
#---

# Combinatorial control of smooth muscle-specific gene expression.
# PMID: 12740224
# Molecular regulation of vascular smooth muscle cell 
# differentiation in development and disease
# PM:15269336
#---


xiaon.counts           <- xiaon.genes %$% 
                          counts %>%
                          DGEList(genes = xiaon.genes$annotation) %>%
                          calcNormFactors() %>% rpkm(log = F)
rownames(xiaon.counts) <- mapIds( org.Mm.eg.db, 
                                      keys      = xiaon.genes$annotation$GeneID %>% 
                                                  as.character(),
                                      column    = 'SYMBOL', 
                                      keytype   = 'ENTREZID', 
                                      multiVals = 'first') %>%
                          make.names(unique = T)
mouse.smc.genes <- list( Smtn = 29856, Acta2 = 11475, Myh11 = 17880,
                         Cnn1 = 12797, Tagln = 21345, Cald1 = 109624)
                         
mouse.smc.counts <- xiaon.counts[names(mouse.smc.genes),1]

#---
# https://www.proteinatlas.org/humanproteome/brain
# 419 brain enriched genes
# Proteins specifically expressed in neurons
#---
human.brain.genes <- list( OPALIN = 93377, GFAP = 2670,OMG = 4974,
                           OLIG2  = 10215, GRIN1 = 2902, NEUROD6 = 63974,
                           SLC17A7  = 57030,CREG2   = 200407, NEUROD2 = 4761,
                           C1orf61  = 10485,ZDHHC22 = 283576, KCNJ9 = 3765)

#---
# no C1orf61 in mouse
#---
mouse.brain.genes <- list( Opalin = 226115, Gfap   = 14580,  Omg = 18377,
                           Olig2  = 50913,  Grin1  = 14810,  Neurod6 = 11922,
                           Slc17a7 = 72961, Creg2  = 263764, Neurod2 = 18013,
                           Zdhhc22 = 238331, Kcnj9 = 16524)

mouse.brain.counts <- xiaon.counts[names(mouse.brain.genes),1]


#---
# TITLE:
# Bioinformatic identification and characterization of 
# human endothelial cell-restricted genes.
# PMID: 20509943
#---

human.endothelium.genes <- list( MMRN1 = 22915,CLDN5 = 7122,VWF   = 7450,
                                 BMX   = 660, ANGPT2 = 285, GJA4  = 2701,
                                 CDH5  = 1003, TIE1  = 7075,ROBO4 = 54538, ECSCR = 641700 )
mouse.endothelium.genes <- list( Mmrn1 = 70945, Cldn5  = 12741, Vwf   = 22371,
                                  Bmx   = 12169, Angpt2 = 11601, Gja4  = 14612,
                                  Cdh5  = 12562, Tie1   = 21846, Robo4 = 74144, Ecscr = 68545)

mouse.endothelium.counts <- xiaon.counts[names(mouse.endothelium.genes),1]


# now the heart specific genes
# imported from my own curation
#---
setwd('E:\\FuWai\\dead.dir\\wangli.lab\\Others')
mouse.heart.df <- read_excel('cardio_fibro_manual_markers.xlsx', sheet = 2)
mouse.heart.counts <- xiaon.counts[mouse.heart.df$symbol,1]
#---
# the RNA-seq count data follows Count-based (negative binomial)
# distribution. I called the function in edgeR package to test
# if the two groups have different distribution parameter.
# exactTest
#--- 

tissue.genes.df <- data.frame( tissue.counts = c( mouse.brain.counts,
                                                  mouse.smc.counts,
                                                  mouse.endothelium.counts,
                                                  mouse.heart.counts),
                               tissue.source = c( rep(1,length(mouse.brain.counts)),
                                                  rep(2,length(mouse.smc.counts)),
                                                  rep(3,length(mouse.endothelium.counts)),
                                                  rep(4,length(mouse.heart.counts)) ) )
                               
                               
tissues.ggplot <- ggplot(data = tissue.genes.df, aes( x = as.factor(tissue.source),
                                                      y = tissue.counts)) +
                  geom_boxplot() +
                  geom_dotplot( aes(as.factor(tissue.source)), 
                                binaxis = 'y', stackdir = 'center',
                                dotsize = .6, position = position_dodge(0.4)) +
                  xlab('tissue specific genes') +
                  ylab('gene normalized count (FPKM)') +
                  scale_x_discrete(labels = c('brain','smooth muscle','endothelium','heart'))
(tissues.ggplot)

# exactTest(mouse.smc.counts,mouse.brain.counts)




#---
# this is the QC step to parse the 
# ribosome percentage in each library
# the inference method is Picard, see
# the scripts:
# quality_control.sh
# the above script get the picard output
# results
# and this script:
# QCparser.pl
# this Perl script output xiaon.ribosome.csv
#---

setwd('E:\\FuWai\\Bioinfo\\chenlab')
xiaon.ribosome.pct <- read_delim( 'xiaon.ribosome.csv', 
                                  col_names = FALSE,
                                  delim = '\t', col_types = 'cn')
xiaon.ribo.ggplot <- ggplot(data = xiaon.ribosome.pct, aes(x = as.factor(X1), y = X2)) +
                     geom_bar(stat = 'identity') +
                     xlab('biological sample names') +
                     ylab('ribosome percentage') +
                     geom_hline(yintercept = 0.01, color = 'green') +
                     theme( axis.text.x = element_text(angle = 60,hjust = 1, size = 6))

(xiaon.ribo.ggplot)
#---
# now complete the pre-QC and post-QC in ning.xiao project
# let start the limma analysis or other DEG analysis
# aim:
# to find the genes which is responsible for conditioned on blood hypertension
#
#---


xiaon.groups              <- factor( c( rep(1,4), rep(2,5), rep(3,4),
                                         rep(4, 4),rep(5, 4),rep(6, 5) ), 
                                     levels = 1:6,
                                     labels = c( 'AL_HU','AL_ICH','AL_LU', 
                                                  'A', 'A_LU','Veh' ))
xiaon.design              <- model.matrix(~ 0 + xiaon.groups)

colnames(xiaon.design)    <- levels(xiaon.groups)
xiaon.contrast.matrix     <- makeContrasts( A_LU - A, AL_LU - AL_ICH,
                                            levels = xiaon.design)
xiaon.limma.result        <- xiaon.genes %$% counts %>% 
                             DGEList(genes = xiaon.genes$annotation) %>% 
                             calcNormFactors() %>% 
                             voom(design = xiaon.design) %>%
                             lmFit(xiaon.design) %>% 
                             contrasts.fit(xiaon.contrast.matrix ) %>%
                             eBayes()
mm10.gene.symbols         <- mapIds( org.Mm.eg.db, 
                                     keys      = xiaon.genes$annotation$GeneID %>% 
                                                 as.character(),
                                     column    = 'SYMBOL', 
                                     keytype   = 'ENTREZID', 
                                     multiVals = 'first') %>%
                             make.names(unique = T)

xiaon.result.1            <- topTable( xiaon.limma.result , 
                                            coef          = 1,
                                            number        = Inf, 
                                            adjust.method = 'BH', 
                                            sort.by       = 'none') %>%
                              cbind(mm10.gene.symbols) %>%
                              arrange(P.Value) %>%
                              dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
                              dplyr::rename(Symbol = mm10.gene.symbols)

xiaon.result.2            <- topTable( xiaon.limma.result , 
                                            coef          = 2,
                                            number        = Inf, 
                                            adjust.method = 'BH', 
                                            sort.by       = 'none') %>%
                              cbind(mm10.gene.symbols) %>%
                              arrange(P.Value) %>%
                              dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
                              dplyr::rename(Symbol = mm10.gene.symbols)
 

#---
# export the result to two sheets
#
#---

setwd('E:\\FuWai\\Bioinfo\\chenlab')
xiaon.wb <- createWorkbook()

addWorksheet(xiaon.wb, 'A_LU - A')
addWorksheet(xiaon.wb, 'AL_LU - AL_ICH')


writeData(xiaon.wb, sheet = 1, xiaon.result.1 )
writeData(xiaon.wb, sheet = 2, xiaon.result.2 )
saveWorkbook(xiaon.wb, 'xiaon.2018-01-17.xlsx', overwrite = TRUE)


#---
# to check the phenotype hypothesis
# if the cell turn into cardiac cell?
# this is random effect
#---

xiaon.1.P.value         <- xiaon.result.1$P.Value
names(xiaon.1.P.value)  <- xiaon.result.1$GeneID
xiaon.1.P.value         <- sort(xiaon.1.P.value, decreasing = TRUE)

heart.core.genes        <- data.frame( diseaseId = rep('cardioCore',length(mouse.heart.df$geneID)), 
                                       geneId    = mouse.heart.df$geneID, check.names = TRUE)

xiaon.1.GSEA              <- GSEA( xiaon.1.P.value, TERM2GENE = heart.core.genes, 
                                   nPerm = 10^5, maxGSSize = 5000, pvalueCutoff = 1)

result.GSEA.table       <- grid.newpage() %>% {summary(core2.GSEA)} %>% as.data.frame %>%
                           dplyr::select(-(c(core_enrichment, leading_edge, ID))) %>%
                           tableGrob(rows = NULL) %>%
                           grid.draw()

gseaplot(xiaon.1.GSEA , geneSetID = 'cardioCore')


#---
# KEGG enrichment analysis
#---

xiaon.kegg.tidy      <- xiaon.result.1 %>% 
                        filter(P.Value < 0.05) %>%
                        dplyr::select(GeneID) %>%
                        unlist() %>%
                        enrichKEGG( ., 
                                    organism = 'mouse', 
                                    pvalueCutoff  = 0.05, 
                                    pAdjustMethod = 'none',
                                    qvalueCutoff  = 1) %>%
                        summary() %$%
                        {data.frame( kegg.pvalue  = -log(pvalue),
                                     kegg.pathway = Description )}
xiaon.kegg.ggplot   <- xiaon.kegg.tidy[1:15,] %>% 
                       ggplot( aes( x = reorder(kegg.pathway, kegg.pvalue), 
                                                   y = kegg.pvalue)) + 
                       geom_bar( stat = 'identity', width = 0.6, 
                                 position = position_dodge(width = 0.05), size = 25) +
                       theme( axis.text.x = element_text(angle = 60,hjust = 1, size = 8),
                              axis.text.y = element_text(hjust = 1, size = 10)) +
                       ylab('-log(pvalue)') + 
                       xlab('xioan.KEGG.plot') + 
                       coord_flip()

xiaon.tableGrob     <- tableGrob( xiaon.kegg.tidy, 
                                  theme = ttheme_default( 
                                            base_size = 7, 
                                            padding   = unit(c(0.5, 0.5), 'mm'),
                                            plot.margin = unit(c(1,0,1,0),'mm')
                                            ) )
ggdraw() + draw_plot(xiaon.tableGrob, x = 0, y = 0)
