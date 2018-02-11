# @author Yisong Zhen
# @since  2017-06-06
# @update 2018-02-09
# @parent 
#    processed results are from the cardiac_cell_ChIP.sh
#    the above shell script have two branches, master and bowpic
#    which use different peak calling program. MACS2 & epic
#
#--- 


#---
# Public data source
# Series GSE52386 
# Nord AS, Blow MJ, Attanasio C, Akiyama JA et al. 
# Rapid and pervasive changes in genome-wide enhancer usage \
# during mammalian development. Cell 2013 Dec 19;155(7):1521-31. \
# PMID: 24360275
#---


#---
# normalization?
# https://www.biostars.org/p/195689/
# https://www.biostars.org/p/42291/
#----

#---
# how to install R lib in the local. I
# I do not try this, instead, I use the miniconda
# to manage the software installation and version control
#
# brainchild from
# http://genomicsclass.github.io/book/pages/read_counting.html
#
# .libPaths('/home/zhenyisong/data/rlib/')
# chooseCRANmirror()
# install.packages( 'tidyverse',lib ='/wa/zhenyisong/rlib')
# install.packages( 'tidyverse', repo = 'http://mirrors.ustc.edu.cn/CRAN/')
# source("http://bioconductor.org/biocLite.R")
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/") 
# biocLite(pkgs, lib = '/home/zhenyisong/data/rlib')
#---


pkgs <- c( 'tidyverse', 'GenomicRanges',
           'ChIPseeker', 'rtracklayer',
           'GenomicAlignments', 'BiocParallel',
           'Rsamtools','magrittr',
           'TxDb.Mmusculus.UCSC.mm10.knownGene',
           'Mus.musculus', 'ggbio',
           'BSgenome.Hsapiens.UCSC.hg38',
           'BSgenome.Hsapiens.UCSC.hg38.Rbowtie',
           'BSgenome.Mmusculus.UCSC.mm10',
           'BSgenome.Mmusculus.UCSC.mm10.Rbowtie',
           'BSgenome.Rnorvegicus.UCSC.rn6',
           'BSgenome.Rnorvegicus.UCSC.rn6.Rbowtie')

load.lib <- lapply(pkgs, require, character.only = TRUE)


working.env <- 'window'
linux.path  <- file.path('/wa/zhenyisong/cardiodata/GSE52386/bwa')
window.path <- file.path('D:/yisong.data')

#---
# [zhenyisong@mgt data]$ md5sum xiaon.Rdata
# 954e0ee9a5969a873f2f0fd9791d0375  xiaon.Rdata
#---

switch( working.env, linux  = {   setwd(linux.path);
                                  if( file.exists('cellChIP.Rdata') ) {
                                      load('cellChIP.Rdata');
                                  }
                              } ,
                     window = {   setwd(window.path);
                                  if( file.exists('cellChIP.Rdata') ) {
                                      load('cellChIP.Rdata');
                                  }
                              }  )


GSE52386.data.path   <- file.path('/wa/zhenyisong/cardiodata/GSE52386/bwa')

#---
# aim: the function to read-in the macs peak calling result
# @param
#       the peaking result in xls format from macs2 output.
# @return
#      the GRange object which store all the data
#---
read.macs2.func      <- . %>% read.delim( header = TRUE, sep = '\t',
                                          fill   = TRUE, comment.char = '#', 
                                          stringsAsFactors = FALSE) %$% 
                              { GRanges(  
                                  seqname      = chr,
                                  ranges       = IRanges( start = start, 
                                                          end   = end),
                                  peak.length  = length,
                                  abs.summit   = abs_summit,
                                  pileup       = pileup,
                                  log.pvalue   = X.log10.pvalue.,
                                  foldChange   = fold_enrichment,
                                  logqvalue    = X.log10.qvalue.,
                                  macs2.name   = name) }

#---
# https://support.bioconductor.org/p/83599/
# to merge the GRange object?
#---

GSE52386.macs2.features           <- list.files( GSE52386.data.path, 
                                                 pattern = '_peaks.xls') %>%
                                     map(read.macs2.func)

GSE52386.macs2.reduced.features   <- GSE52386.macs2.features %>% 
                                     do.call(getMethod(c, "GenomicRanges"), .) %>%
                                     IRanges::reduce()


#---
# yieldSize should have big value, otherwise the system will throw out 
# a exception :
# Error: 'bplapply' receive data failed:
# Error: unexpected symbol in "Error: 'bplapply' receive"
#---

GSE52386.bam.filenames     <- list.files(GSE52386.data.path, pattern = '*.bam$') %>% 
                            {.[seq(1,length(.), 2)]} %>% .[c(1,13:18)]
#setwd(GSE52386.data.path)
indexBam(GSE52386.bam.filenames)
GSE52386.bamFileList       <- Rsamtools::BamFileList(
                                         GSE52386.bam.filenames, 
                                         yieldSize = 7500000)

                            

#names(GSE52386.macs2.bams)  <- paste0('macs2_result','_',1:7,sep = '')
names(GSE52386.macs2.bams)  <- c( 'Heart_E11.5','Heart_E14.5','Heart_E17.5',
                                  'Heart_P0','Heart_P7','Heart_P21','Heart_P56')

#---
# https://support.bioconductor.org/p/84541/
# https://github.com/genomicsclass/labs/blob/master/
# biocadv_6x/bioc2_parallel.Rmd#implicit-parallelism-through-biocparallel
# BiocParallel and implicit parallelization
# Concurrent counting of RNA-seq reads
# labs/biocadv_6x/bioc2_parallel.Rmd
#---

register( MulticoreParam( workers = 4) )
GSE52386.read.counts <- summarizeOverlaps(  
                           features      = GSE52386.macs2.reduced.features,
                           reads         = GSE52386.macs2.bams, 
                           ignore.strand = TRUE, 
                           singleEnd     = TRUE ) %>% assay()



# visiolization
#---
mouse.txdb    <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(Mus.musculus)
gene.positive.control <- genes( mouse.txdb, 
                                filter = list(gene_id = 16001)) %>%
                         keepSeqlevels('chr7')

pos.control.reads   <- . %>% 
                       readGAlignments( param = ScanBamParam( 
                                                   which = gene.positive.control), 
                                        use.names = T)
GSE52386.bam.names  <- list.files(GSE52386.data.path, pattern = '*.bam$') %>% 
                       {.[seq(1,length(.), 2)]} %>% .[c(1,13:18)]

GSE52386.readBam.results <- map(GSE52386.bam.names, pos.control.reads)


gene.model   <- ggbio::autoplot( Mus.musculus, which = gene.positive.control, 
                          columns = c('GENENAME', 'SYMBOL'), 
                          names.expr = 'GENENAME::SYMBOL')
thocs5.cov   <- ggbio::autoplot( GSE52386.readBam.results[[1]], 
                                 which = gene.positive.control) + 
                       ylim(0,100) + ylab('thocs5')
input.cov    <- ggbio::autoplot( GSE52386.readBam.results[[2]], 
                                 which = gene.positive.control ) + 
                       ylim(0,100) + ylab('input')
ggbio::tracks(gene.model, thocs5.cov, input.cov, heights = c(1,2,2))




setwd(GSE52386.data.path)
save.image('cellChIP.Rdata')
q('no')
