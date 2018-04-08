# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-04-08
# @parent  blood_pipeline_QC.sh
#---

pkgs                   <- c( 'tidyverse', 'GenomicRanges',
                             'ChIPseeker', 'rtracklayer',
                             'GenomicAlignments', 'BiocParallel',
                             'Rsamtools','magrittr', 'DESeq2',
                             'stringr', 'JASPAR2016', 'TFBSTools',
                             'seqLogo', 'RSQLite', 'DBI',
                             'TxDb.Dmelanogaster.UCSC.dm6.ensGene',
                             'ggbio', 'ChIPpeakAnno',
                             'BSgenome.Dmelanogaster.UCSC.dm6')
         
load.lib               <- lapply(pkgs, require, character.only = TRUE)
  
macs2.ChIRP.path       <- file.path('/wa/zhenyisong/results/chenlab/songli/pipelineQC/bwa')

setwd(macs2.ChIRP.path)
rox2.pearson.results   <- read.delim( 'merge_peaks.xls.corr.xls',
                                       header = TRUE, sep = '\t',
                                       fill   = TRUE, comment.char = '#', 
                                       stringsAsFactors = FALSE) %>%
                          filter(fold_enrichment >= 2) %>%
                          filter(correlation >= 0.3)   %>%
                          filter(aver_coverage >= 1.5)
rox2.pearson.GR        <- rox2.pearson.results %$%
                         { GRanges(  
                             seqname  = as.character(chr),
                             ranges   = IRanges( start = start, 
                                                 end   = end ) ) } %>%
                         .[width(.) >= 2300]

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
                            macs2.name   = name ) }
macs2.ChIRP.features   <- list.files( macs2.ChIRP.path, 
                                      full.names = TRUE,
                                      pattern    = '_peaks.xls') %>%
                          map(read.macs2.func)

even.ChIRP.dm6.macs2  <- macs2.ChIRP.features[[1]]
odd.ChIRP.dm6.macs2   <- macs2.ChIRP.features[[2]]
roX2.full.macs2       <- macs2.ChIRP.features[[3]]
ChIRP.dm6.annot       <- toGRanges( TxDb.Dmelanogaster.UCSC.dm6.ensGene, feature = 'gene')
## keep the seqnames in the same style
seqlevelsStyle(roX2.full.macs2)        <- seqlevelsStyle(ChIRP.dm6.annot)
## do annotation by nearest TSS       
roX2.peaks.annot                       <- annotatePeakInBatch( roX2.full.macs2,
                                                AnnotationData = ChIRP.dm6.annot)
summary(mcols(roX2.peaks.annot)$foldChange)
sum( mcols(roX2.peaks.annot)$peak.length > 2300 & 
     mcols(roX2.peaks.annot)$foldChange > 2)
rox2.filter <- mcols(roX2.peaks.annot)$peak.length > 2300 & 
               mcols(roX2.peaks.annot)$foldChange > 2
rox2.sex.chromosome <- seqnames(roX2.peaks.annot) == 'chrX'
roX2.peaks.annot[rox2.sex.chromosome]

seqlevelsStyle(even.ChIRP.dm6.macs2)   <- seqlevelsStyle(ChIRP.dm6.annot)
even.peaks.annot                       <- annotatePeakInBatch( even.ChIRP.dm6.macs2,
                                                AnnotationData = ChIRP.dm6.annot)
seqlevelsStyle(odd.ChIRP.dm6.macs2)    <- seqlevelsStyle(ChIRP.dm6.annot)
odd.peaks.annot                        <- annotatePeakInBatch( odd.ChIRP.dm6.macs2,
                                                AnnotationData = ChIRP.dm6.annot)