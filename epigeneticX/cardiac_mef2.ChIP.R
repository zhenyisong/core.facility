#--
# @author Yisong Zhen
# @since  2018-09-18
# @update 2018-09-20
# @parent cardiac_mef2.ChIP.py
#---

#---
# 
# @raw_data_path
# rat, SE
# two inputs, two controls, = one replicate
#---

pkgs     <- c( 'tidyverse','TFBSTools', 'PWMEnrich',
               'BSgenome.Rnorvegicus.UCSC.rn6','org.Rn.eg.db',
               'magrittr')
load.lib <- lapply(pkgs, require, character.only = TRUE)

#---
# define constants
#---
raw_data_path <- file.path('/wa/zhenyisong/cardiodata/PRJEB23434')
rn6.genome    <- BSgenome.Rnorvegicus.UCSC.rn6


#---
# define functions
#---

.get_macs2_result_filenames <- function( datadir, file_pattern = '.xls') {
      macs2_filenames <- list.files( path = datadir, 
                                     full.names = TRUE,
                                     pattern = file_pattern)
}

.get_macs_peaks_result <- function(macs.filename) {
    macs2.result     <- read.delim( macs.filename, 
                                    comment.char = '#') %$%
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
    return(macs2.result)
}

.get_macs2_peak_seqs <- function(seqs_granges, genome = rn6.genome) {
   macs2_seqs <- getSeq(genome, seqs_granges)
   return(macs2_seqs)
}

macs2.rat.seqs <-  .get_macs2_result_filenames(raw_data_path) %>%
                   map(.get_macs_peaks_result) %>%
                   map(.get_macs2_peak_seqs)
mef2.enriched  <- motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)