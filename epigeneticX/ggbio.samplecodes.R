#---
# Introduction to ggbio
# 
# Leonardo Collado-Torres
# October 8, 2013
# 
# ggbio
# :  visualization toolkits for genomic data
#
# ggbio - Visualize genomic data 
# http://www.sthda.com/english/wiki/ggbio-visualize-genomic-data
#---
library('ggbio')
load(system.file( 'data', 'hg19IdeogramCyto.rda', 
                   package  = 'biovizBase', 
                   mustWork = TRUE))
p.ideo <- plotIdeogram(hg19IdeogramCyto, 'chr1')
print(p.ideo)

library('GenomicRanges')
data(hg19Ideogram, package= 'biovizBase')
autoplot(seqinfo(hg19Ideogram)[paste0('chr', 1:13)])

## Generate
gr <- GRanges( sample(paste0('chr', 1:13), 50, TRUE), 
               ranges = IRanges(round(runif(50, 1, 1e8)), 
               width  = 1000))
## Get chr lengths
seqlengths(gr) <- seqlengths(hg19Ideogram)[names(seqlengths(gr))]
## Add a group variable
gr$group <- factor(sample(letters[1:4], 50, TRUE))

autoplot(seqinfo(gr)) + layout_karyogram( gr, aes( fill = group, color = group))


y  <- as.vector(sapply(1:5, function(z) rnorm(100, 10 * z * (-1)^z, 2 * z)))
df <- data.frame(y = y + 1.1 * abs(min(y)), x = seq_len(length(y)))

p.noisy <- ggplot(df, aes(x=x, y=y)) + geom_line()
p.noisy

exons  <- GRanges(rep('chr1', 2), IRanges(c(101, 301), width = 100))
p.exon <- autoplot(exons)
print(p.exon)


final <- tracks( p.ideo, 'Coverage' = p.noisy, 
                         'Exons'    = p.exon, 
                          heights   = c(2, 6, 3), 
                          title     = 'Tracks plot!') + 
         ylab('') + 
         theme_tracks_sunset()
print(final)


library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggbio)
library(biovizBase)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
columns <- c('CHR','CHRLOC','CHRLOCEND')
sel     <- select( org.Hs.eg.db, 'ENSG00000119541', columns,
                   keytype = 'ENSEMBL')
wh <- GRanges('chr18', IRanges(63389191,63422519), strand = Rle('-'))
p1 <- autoplot(txdb, which = wh, names.expr = 'gene_id')
p2 <- autoplot(txdb, which = wh, stat = 'reduce', color = 'brown', fill = 'brown')
p.ideo <- getIdeogram('hg38', cytoband = TRUE, subchr = 'chr18')
tracks(p.ideo, fill = p1, reduce = p2, heights = c(1.5,5,1)) + 
             ylab('')  + 
             xlim(wh)  +
             theme_tracks_sunset()