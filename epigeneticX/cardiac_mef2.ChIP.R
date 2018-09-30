#--
# @author Yisong Zhen
# @since  2018-09-18
# @update 2018-09-30
# @parent cardiac_mef2.ChIP.py
#---

#---
# 
# @raw_data_path
# rat, SE
# two inputs, two controls, = one replicate
#---

pkgs     <- c( 'tidyverse','TFBSTools', 'PWMEnrich','rtracklayer',
               'BSgenome.Rnorvegicus.UCSC.rn6','org.Rn.eg.db',
               'magrittr','motifcounter','MotifDb','ggseqlogo',
               'BCRANK','DiffLogo')
load.lib <- lapply(pkgs, require, character.only = TRUE)

#---
# define constants
#---
raw_data_path <- file.path('/wa/zhenyisong/cardiodata/PRJEB23434')
rn6.genome    <- BSgenome.Rnorvegicus.UCSC.rn6


#---
# define functions
#---

.get_macs2_result_filenames <- function( datadir, file_pattern = '.bed') {
      macs2.summit.names <- list.files( path = datadir, 
                                        full.names = TRUE,
                                        pattern = file_pattern)
}

.get_macs_summit_result <- function(macs.filename) {
    macs2.summit    <- import(macs.filename, format = 'BED')
    return(macs2.summit)
}

.get_macs2_peak_seqs <- function(seqs_granges, genome = rn6.genome) {
   macs2_seqs <- getSeq(genome, seqs_granges)
   return(macs2_seqs)
}


macs2.rat.seqs <-  .get_macs2_result_filenames(raw_data_path) %>%
                   map(.get_macs_summit_result) %>% 
                   map(function(x) x + 200) %>% 
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
                                   'jolma2013'))[[1]]

motif.mef2.normed   <- as.list( query(query( query( 
                                          MotifDb, 
                                          'hsapiens'), 
                                          'mef2'), 
                                   'jolma2013'))[[1]] %>%
                       normalizeMotif()
ggplot() + geom_logo( motif.mef2 ) + theme_logo()
#seqLogo(motif.mef2)
getwd()
setwd(raw_data_path)
peak.filename <- 'temp.fa'
writeXStringSet(macs2.rat.seqs[[1]], peak.filename , format = 'fasta')
bcrank.out <- bcrank( peak.filename, 
                      restarts = 25, 
                      use.P1   = TRUE, 
                      use.P2   = TRUE)
file.remove(peak.filename)
toptable(bcrank.out)
top.motif     <- toptable(bcrank.out , 1)
weight.matrix <- pwm(top.motif, normalize = FALSE)
ggplot() + geom_logo( motif.mef2 ) + theme_logo()
weight.matrix.normalized <- pwm(top.motif, normalize = TRUE)
seqLogo(weight.matrix.normalized)
#diffLogoFromPwm(pwm1 = motif.mef2 , pwm2 = weight.matrix.normalized)
PWMSimilarity(motif.mef2, weight.matrix.normalized, method = 'Pearson' )
#PFMSimilarity(motif.mef2,weight.matrix)
# save.image('temp.Rdata')
result <- motifEnrichment(macs2.rat.seqs[[1]], motif.mef2, bg.2)