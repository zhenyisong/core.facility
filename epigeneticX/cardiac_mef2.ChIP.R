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
               'magrittr','motifcounter','MotifDb','seqLogo')
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
#mef2.enriched  <- motifEnrichment(sequences, PWMLogn.dm3.MotifDb.Dmel)

#---
# Ideally, the DNA sequence for estimating the background model 
# should be representative (or even the same) as the sequences 
# that are latter analysed (e.g. for motif hit enrichment).
#---
bg.2         <- readBackground(macs2.rat.seqs[[1]], 2)
motif.mef2   <- as.list( query(query( query( 
                                          MotifDb, 
                                          'hsapiens'), 
                                          'mef2'), 
                                   'jolma2013'))[[1]] %>%
                normalizeMotif()
seqLogo(motif.mef2)
result <- motifEnrichment(macs2.rat.seqs[[1]], motif.mef2, bg.2)