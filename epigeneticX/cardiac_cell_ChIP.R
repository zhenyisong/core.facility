# @author Yisong Zhen
# @since  2017-06-06
# @update 2018-02-05
# @parent 
#    processed results are from the cardiac_cell_ChIP.sh
#    the above shell script have two branches, master and bowpic
#    which use different peak calling program. MACS2 & epic
#
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
           'Mus.musculus', 
           'BSgenome.Hsapiens.UCSC.hg38',
           'BSgenome.Hsapiens.UCSC.hg38.Rbowtie',
           'BSgenome.Mmusculus.UCSC.mm10',
           'BSgenome.Mmusculus.UCSC.mm10.Rbowtie',
           'BSgenome.Rnorvegicus.UCSC.rn6',
           'BSgenome.Rnorvegicus.UCSC.rn6.Rbowtie')

load.lib <- lapply(pkgs, require, character.only = TRUE)


working.env <- 'linux'
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
GSE52386.macs2.features   <- list.files( GSE52386.data.path, 
                                        pattern = '_peaks.xls')[1] %>%
                            map(read.macs2.func) %>% 
                            do.call(getMethod(c, "GenomicRanges"), .) %>%
                            IRanges::reduce()



GSE52386.macs2.bams      <- list.files(GSE52386.data.path, pattern = '*.bam$') %>% 
                            {.[seq(1,length(.), 2)]} %>% .[c(1,13:18)] %>%
                            Rsamtools::BamFileList(.)

names(GSE52386.macs2.bams)  <- paste0('macs2_result','_',1:7,sep = '')


# https://support.bioconductor.org/p/84541/
GSE52386.overlap <- summarizeOverlaps(  
                      features      = GSE52386.macs2.features,
                      reads         = GSE52386.macs2.bams, 
                      ignore.strand = TRUE, 
                      singleEnd     = FALSE, 
                      fragments     = TRUE,
                      yieldSize     = 7500000 ) %>% assay()

setwd(GSE52386.data.path)
save.image('cellChIP.Rdata')
q('no')
