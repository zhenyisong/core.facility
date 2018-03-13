# @author Yisong Zhen
# @since  2017-06-06
# @update 2018-03-09
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
           'Rsamtools','magrittr', 'DESeq2',
           'stringr', 'JASPAR2016', 'TFBSTools',
           'seqLogo', 'RSQLite', 'DBI',
           'TxDb.Mmusculus.UCSC.mm10.knownGene',
           'Mus.musculus', 'ggbio', 'ChIPpeakAnno',
           'BSgenome.Hsapiens.UCSC.hg38',
           'BSgenome.Hsapiens.UCSC.hg38.Rbowtie',
           'BSgenome.Mmusculus.UCSC.mm10',
           'BSgenome.Mmusculus.UCSC.mm10.Rbowtie',
           'BSgenome.Rnorvegicus.UCSC.rn6',
           'BSgenome.Rnorvegicus.UCSC.rn6.Rbowtie')

load.lib    <- lapply(pkgs, require, character.only = TRUE)


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
GSE52386.macs2.bams        <- Rsamtools::BamFileList(
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

heart.dev.groups     <- c( 'E11_5','E14_5','E17_5',
                           'P0','P7','P21','P56') %>% 
                        data.frame(heartDevs = as.factor(.))
rownames(heart.dev.groups) <- colnames(GSE52386.read.counts)
heart.ChIP.deseq     <- DESeqDataSetFromMatrix( countData = GSE52386.read.counts,
                                                colData   = heart.dev.groups,
                                                design    = ~ heartDevs)

heart.ChIP.diffbind  <- DESeq(heart.ChIP.deseq)


ChIP.mm10.annot <- toGRanges(TxDb.Mmusculus.UCSC.mm10.knownGene, feature = 'gene')
## keep the seqnames in the same style
seqlevelsStyle(GSE52386.macs2.reduced.features) <- seqlevelsStyle(ChIP.mm10.annot)
## do annotation by nearest TSS
peaks.annot <- annotatePeakInBatch( GSE52386.macs2.reduced.features,
                                    AnnotationData = ChIP.mm10.annot)
#Bmp10 bone morphogenetic protein 10 [ Mus musculus (house mouse) ]
#Gene ID: 12154, updated on 13-Feb-2018
# chr6
BMP10.Granges <- peaks.annot[mcols(peaks.annot)$feature == as.character(12154)]

full.comparisions    <- resultsNames(heart.ChIP.diffbind)
# result <- results(heart.ChIP.diffbind, name = 'heartDevs_P0_vs_E11_5')
targets.int <- names(BMP10.Granges) %>% 
               str_extract( '\\d+') %>% as.integer()
bmp10.inx   <- c()
for(i in 1:length(full.comparisions)) {
    result <- results(heart.ChIP.diffbind, name = full.comparisions[i])
    result <- result %>% as.data.frame() %>%
              mutate(row.inx = 1:nrow(.) )%>%
              .[targets.int,] %>% 
              dplyr::filter(abs(log2FoldChange) >= 3) %>%
              dplyr::select(row.inx) %>% unlist()
    bmp10.inx <- c(bmp10.inx, result)
        
}

bmp10.inx <- unique(bmp10.inx)


#---
# Can I use DESeq2 to analyze a dataset without replicates?
# see DESeq2 working manual
#---



#---
# peak visiolization and manual QC check
# EntrezGeneID: 16001
# Gene Symbol : Igf1r
#---
mouse.txdb    <- TxDb.Mmusculus.UCSC.mm10.knownGene
columns(Mus.musculus)
Igf1r.positive.control <- genes( mouse.txdb, 
                                 filter = list(gene_id = 12154)) %>%
                          keepSeqlevels('chr6')
Igf1r.enhancer.control <- promoters(Igf1r.positive.control, upstream = 80000) 

# class(Igf1r.positive.control)
# > class(Igf1r.positive.control)
# [1] "GRanges"
# attr(,"package")
# [1] "GenomicRanges"
#---
#Igf1r.control.reads   <- . %>% 
#                         readGAlignments( param = ScanBamParam( 
#                                                     which = Igf1r.enhancer.control), 
#                                          use.names = T)
#---
GSE52386.bam.names  <- list.files(GSE52386.data.path, pattern = '*.bam$') %>% 
                       {.[seq(1,length(.), 2)]} %>% .[c(1,13:18)]
#---
# deprecated!!
#GSE52386.readBam.results <- map(GSE52386.bam.names, Igf1r.control.reads)
#gene.model   <- autoplot( Mus.musculus, which = Igf1r.positive.control, 
#                                 columns = c('GENENAME', 'SYMBOL'), 
#                                 names.expr = 'GENENAME::SYMBOL')
#thocs5.cov   <- autoplot( GSE52386.bam.names[1], 
#                          method = 'estimate' ,
#                          aes(y = log(..coverage..)),
#                          which  =  Igf1r.enhancer.control )  
#---

heart.dict          <- list( 'SRR1029001.bam' = 'E11.5', 'SRR1029878.bam' = 'E14.5',
                             'SRR1029880.bam' = 'E17.5', 'SRR1029882.bam' = 'P0',
                             'SRR1029884.bam' = 'P7',    'SRR1029886.bam' = 'P21',
                             'SRR1029888.bam' = 'P56')

heart.coverage.func <- . %>% 
                      autoplot( geom   = 'polygon',
                                size   = 0.5,
                                ylab   = heart.dict[[.]],
                                which  = bmp10.enhancer.grange ) 
heart.coverage.list <- map(GSE52386.bam.names, heart.coverage.func)
tracks( heart.coverage.list)  + ylim(c(1,16))


bmp10.enhancer.grange <- GRanges( seqname = 'chr6', 
                                  ranges  = IRanges( start = 87399000 , 
                                                     end   = 87401000),
                                  strand  = '+') 
mm10.ucsc <- BSgenome.Mmusculus.UCSC.mm10
seqlevelsStyle(mm10.ucsc)
seqlevelsStyle(bmp10.enhancer.grange)
bmp10.enhancer <- BSgenome::getSeq(mm10.ucsc, bmp10.enhancer.grange)


# connect to the sqlite file
#
#jaspar.sqlite     <- dbConnect( drv = SQLite(), 
#                                dbname = JASPAR2016@db)
# get a list of all tables
#jaspar.all.tables <- dbListTables(jaspar.sqlite)
opts  <- list()
opts[['species']] <- 9606
opts[['name']]    <- 'SRF'
PFMatrixList <- getMatrixSet(JASPAR2016, opts)
pwm          <- getMatrixSet(JASPAR2016, opts)[[1]] %>% toPWM()

siteset <- searchSeq( pwm, bmp10.enhancer, 
                      seqname   = 'seq1', 
                      min.score = '60%', 
                      strand    = '*')

#---
# this is a beta module to draw the seqlogo
#---
tbx6 <- rbind(  A = c( 14, 47,  3,  5, 68,  1, 93,   0,  1, 10),
                C = c( 16, 12,  6, 90,  2, 98,  2, 100, 84, 26),
                G = c( 47,  3,  0,  4, 25,  0,  3,   0,  2, 5),
                T = c( 23, 38, 91,  1,  5,  1,  2,   0, 13, 59)) %>%
        {./100}   %>%
        makePWM() %>% 
        seqLogo()

setwd(GSE52386.data.path)
save.image('cellChIP.Rdata')
q('no')
