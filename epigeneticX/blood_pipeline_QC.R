# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-04-10
# @parent  blood_pipeline_QC.sh
#---

pkgs                           <- c( 'tidyverse', 'GenomicRanges',
                                     'ChIPseeker', 'rtracklayer',
                                     'GenomicAlignments', 'BiocParallel',
                                     'Rsamtools','magrittr', 'DESeq2',
                                     'stringr', 'JASPAR2016', 'TFBSTools',
                                     'seqLogo', 'RSQLite', 'DBI',
                                     'TxDb.Dmelanogaster.UCSC.dm6.ensGene',
                                     'ggbio', 'ChIPpeakAnno',
                                     'BSgenome.Dmelanogaster.UCSC.dm6')
         
load.lib                       <- lapply(pkgs, require, character.only = TRUE)
  
replication.rox2.ChIRP.path    <- file.path(
                                  '/wa/zhenyisong/results/chenlab/songli/pipelineQC/bowtie1')

setwd(replication.rox2.ChIRP.path)
#---
# this is the replication result for the 
# Mol. Cell data and their result.
# using bowtie1 and macs14
#---
rox2.defacto.pearson.results   <- read.delim( 'merge_peaks.xls.corr.xls',
                                               header = TRUE, sep = '\t',
                                               fill   = TRUE, comment.char = '#', 
                                               stringsAsFactors = FALSE) %>%
                                  filter(fold_enrichment >= 2) %>%
                                  filter(correlation >= 0.3)   %>%
                                  filter(aver_coverage >= 1.5) %$%
                                  { GRanges(  
                                      seqname  = as.character(chr),
                                      ranges   = IRanges( start = start, 
                                                          end   = end ) ) } %>%
                                  .[width(.) >= 2300]

#---
# this result is the bwa program with
# my own understanding, and own protocol
# using macs2
# this method does not replicate the 
# true result and I therefore check the overlapp
#---
rox2.change.ChIRP.path    <- file.path(
                                  '/wa/zhenyisong/results/chenlab/songli/pipelineQC/bwa')

setwd(rox2.change.ChIRP.path)

rox2.private.pearson.results   <- read.delim( 'rox2_peaks.xls.corr.xls',
                                               header = TRUE, sep = '\t',
                                               fill   = TRUE, comment.char = '#', 
                                               stringsAsFactors = FALSE) %>%
                                  filter(fold_enrichment >= 2) %>%
                                  filter(correlation >= 0.3)   %>%
                                  filter(aver_coverage >= 1.5) %$%
                                  { GRanges(  
                                      seqname  = as.character(chr),
                                      ranges   = IRanges( start = start, 
                                                          end   = end ) ) } %>%
                                  .[width(.) >= 2300]

rox2.macs2.xls.results   <- read.delim( file = 'rox2_peaks.xls', 
                                        header = TRUE, sep = '\t',
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
# the results of this step confirms that 
# the two protocol, macs2 + bwa or macs14 + bowtie1 are
# comparable and interchangble.
#---

two.studies <- suppressWarnings( findOverlaps( query   = rox2.defacto.pearson.results, 
                                               subject = rox2.macs2.xls.results, 
                                               type    = 'any') )


rox2.minor.pearson.results   <- read.delim( file   = 'merge_z_peaks.xls.corr.xls',
                                            header = TRUE, sep = '\t',
                                            fill   = TRUE, comment.char = '#', 
                                            stringsAsFactors = FALSE) %>%
                                filter(fold_enrichment >= 2) %>%
                                filter(correlation >= 0.3)   %>%
                                filter(aver_coverage >= 1.5) %$%
                                { GRanges(  
                                    seqname  = as.character(chr),
                                    ranges   = IRanges( start = start, 
                                                        end   = end ) ) } %>%
                                .[width(.) >= 2300]

