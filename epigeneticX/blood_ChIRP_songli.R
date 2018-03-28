# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-03-28
# @parent  blood_ChIRP_songli.sh
#---

pkgs        <- c( 'tidyverse', 'GenomicRanges',
                  'ChIPseeker', 'rtracklayer',
                  'GenomicAlignments', 'BiocParallel',
                  'Rsamtools','magrittr', 'DESeq2',
                  'stringr', 'JASPAR2016', 'TFBSTools',
                  'seqLogo', 'RSQLite', 'DBI',
                  'TxDb.Hsapiens.UCSC.hg38.knownGene',
                  'Homo.sapiens', 'ggbio', 'ChIPpeakAnno',
                  'BSgenome.Hsapiens.UCSC.hg38',
                  'BSgenome.Hsapiens.UCSC.hg38.Rbowtie')
       
load.lib    <- lapply(pkgs, require, character.only = TRUE)

macs2.ChIRP.path     <- file.path('/wa/zhenyisong/results/chenlab/songli/bwa')

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
#---
# https://support.bioconductor.org/p/83599/
# to merge the GRange object?
#---

macs2.ChIRP.features   <- list.files( macs2.ChIRP.path, 
                                      pattern = '_peaks.xls') %>%
                          map(read.macs2.func)
even.ChIRP.hg38.macs2  <- macs2.ChIRP.features[[1]]
odd.ChIRP.hg38.macs2   <- macs2.ChIRP.features[[2]]
ChIRP.hg38.annot       <- toGRanges( TxDb.Hsapiens.UCSC.hg38.knownGene, feature = 'gene')
## keep the seqnames in the same style
seqlevelsStyle(even.ChIRP.hg38.macs2) <- seqlevelsStyle(ChIRP.hg38.annot)
## do annotation by nearest TSS
even.peaks.annot <- annotatePeakInBatch( even.ChIRP.hg38.macs2,
                                         AnnotationData = ChIRP.hg38.annot)

odd.peaks.annot  <- annotatePeakInBatch( odd.ChIRP.hg38.macs2,
                                         AnnotationData = ChIRP.hg38.annot)