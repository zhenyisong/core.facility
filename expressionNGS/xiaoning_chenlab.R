# this will be the working code to process the 
# RNA-seq data from Xiao'ning who is from Jingzhou lab
# @author Yisong
# @since  2017-12-27
# @update 2018-01-04 
#---

#---
# https://stackoverflow.com/questions/18306362/run-r-script-from-command-line
#---
pkgs <- c( 'tidyverse','Rsubread','org.Mm.eg.db','edgeR',
           'limma', 'DESeq2', 'genefilter','grid',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'cowplot',
           'BSgenome.Mmusculus.UCSC.mm10',
           'BSgenome.Mmusculus.UCSC.mm10.Rbowtie')

load.lib <- lapply(pkgs, require, character.only = TRUE)

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
xioan.sample.names    <- basename(read.1.files) %>% 
                         sub(pattern = '_R1.fq.gz', replacement = '')



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



