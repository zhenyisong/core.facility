# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-04-03
# @parent  blood_ChIRP_songli.sh
#---

pkgs              <- c( 'tidyverse', 'GenomicRanges',
                        'ChIPseeker', 'rtracklayer',
                        'GenomicAlignments', 'BiocParallel',
                        'Rsamtools','magrittr', 'DESeq2',
                        'stringr', 'JASPAR2016', 'TFBSTools',
                        'seqLogo', 'RSQLite', 'DBI',
                        'TxDb.Hsapiens.UCSC.hg38.knownGene',
                        'Homo.sapiens', 'ggbio', 'ChIPpeakAnno',
                        'org.Hs.eg.db',
                        'BSgenome.Hsapiens.UCSC.hg38',
                        'BSgenome.Hsapiens.UCSC.hg38.Rbowtie')
             
load.lib          <- lapply(pkgs, require, character.only = TRUE)

macs2.ChIRP.path  <- file.path('/wa/zhenyisong/results/chenlab/songli/bwa')

read.macs2.func   <- . %>% read.delim( header = TRUE, sep = '\t',
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
                                      pattern = '_peaks.xls')[c(2,4,7)] %>%
                          map(read.macs2.func)
setwd(macs2.ChIRP.path)
even.ChIRP.hg38.macs2  <- macs2.ChIRP.features[[1]]
odd.ChIRP.hg38.macs2   <- macs2.ChIRP.features[[2]]
blood.ChIRP.hg38.macs2 <- macs2.ChIRP.features[[3]]
ChIRP.hg38.annot       <- toGRanges( TxDb.Hsapiens.UCSC.hg38.knownGene, feature = 'gene')

seqlevelsStyle(blood.ChIRP.hg38.macs2) <- seqlevelsStyle(ChIRP.hg38.annot)
## do annotation by nearest TSS
blood.peaks.annot <- annotatePeakInBatch( blood.ChIRP.hg38.macs2,
                                          AnnotationData = ChIRP.hg38.annot)

## keep the seqnames in the same style
seqlevelsStyle(even.ChIRP.hg38.macs2) <- seqlevelsStyle(ChIRP.hg38.annot)
## do annotation by nearest TSS
even.peaks.annot <- annotatePeakInBatch( even.ChIRP.hg38.macs2,
                                         AnnotationData = ChIRP.hg38.annot)

odd.peaks.annot  <- annotatePeakInBatch( odd.ChIRP.hg38.macs2,
                                         AnnotationData = ChIRP.hg38.annot)

# read the nova bam pearson analysis result
# see the blood_ChIRP_songli.sh
# the filtering critera were following 
# the recommendation by Mol. Cells
#
#---
setwd(macs2.ChIRP.path)
nova.pearson.results <- read_tsv('blood_nova_peaks.xls.corr.xls') %>%
                        filter(fold_enrichment >= 2) %>%
                        filter(correlation >= 0.3)   %>%
                        filter(aver_coverage >= 1.5)

z.pearson.results    <- read.delim( 'blood_peaks.xls.corr.xls',
                                    header = TRUE, sep = '\t',
                                    fill   = TRUE, comment.char = '#', 
                                    stringsAsFactors = FALSE) %>%
                        filter(fold_enrichment >= 2) %>%
                        filter(correlation >= 0.3)   %>%
                        filter(aver_coverage >= 1.5)
chromosome.set       <- paste('chr',c(1:22,'X','Y'), sep = '')  

nova.pearson.GR      <- nova.pearson.results %$%
                        { GRanges(  
                            seqname  = paste('chr',chr, sep = '') %>%
                                       as.character(),
                            ranges   = IRanges( start = start, 
                                                end   = end ) ) } %>%
                        {.[seqnames(.) %in% chromosome.set]}
z.pearson.GR         <- z.pearson.results %$%
                        { GRanges(  
                            seqname  = as.character(chr),
                            ranges   = IRanges( start = start, 
                                                end   = end ) ) } %>%
                        {.[seqnames(.) %in% chromosome.set]}


two.studies <- suppressWarnings( findOverlaps( query   = nova.pearson.GR, 
                                               subject = z.pearson.GR, 
                                               type    = 'any') )
               
peaks.after.filter.annot <- annotatePeakInBatch( z.pearson.GR,
                                    AnnotationData = ChIRP.hg38.annot)
peaks.hg38.gene.symbols  <- mapIds( org.Hs.eg.db, 
                                     keys      = peaks.after.filter.annot %>%
                                                 mcols() %$% feature %>%
                                                 unique(),
                                     column    = 'SYMBOL', 
                                     keytype   = 'ENTREZID', 
                                     multiVals = 'first') 
mcols(peaks.after.filter.annot)$feature %>% unique