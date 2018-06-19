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