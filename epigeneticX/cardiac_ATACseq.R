# @author  Yisong Zhen
# @since   2018-06-11
# @update  2018-06-11
# @parent  cardiac_ATACseq.sh
#---

pkgs              <- c( 'tidyverse', 'GenomicRanges',
                        'GenomicAlignments', 'BiocParallel',
                        'Rsamtools','magrittr', 'DESeq2',
                        'stringr', 'cluster','factoextra',
                        'TxDb.Mmusculus.UCSC.mm10.knownGene',
                        'Mus.musculus', 'ggbio', 'ChIPpeakAnno',
                        'org.Mm.eg.db', 'clusterProfiler',
                        'edgeR','DESeq2', 'openxlsx', 'pheatmap','RColorBrewer',
                        'BSgenome.Mmusculus.UCSC.mm10')
             
load.lib          <- lapply(pkgs, require, character.only = TRUE)

raw.data.paths    <- file.path('/wa/zhenyisong/cardiodata/SRP101479/bwa')
ATAC.peak.files   <- list.files( path         = raw.data.paths, 
                                 pattern      = '_peaks.xls$', 
                                 all.files    = FALSE, 
                                 full.names   = TRUE, 
                                 recursive    = FALSE, 
                                 ignore.case  = FALSE, 
                                 include.dirs = FALSE)

read_peaks_xls_file <- function(filename) {
    peak_file       <- read.delim( filename, 
                                   comment.char = '#')
}

change_df_2_grange  <- function(data_df) {
    data_df %$% { GRanges(  
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
}

get_grange_from_df <- function(data_df) {
    data_df %>% read_peaks_xls_file() %>% change_df_2_grange()
}

complete_ATAC_Grange.data <- map(ATAC.peak.files, get_grange_from_df)

common_Granges <- Reduce(subsetByOverlaps,complete_ATAC_Grange.data)

ATAC.peak.bams   <- list.files( path          = raw.data.paths, 
                                 pattern      = '.dedup.sorted.bam$', 
                                 all.files    = FALSE, 
                                 full.names   = TRUE, 
                                 recursive    = FALSE, 
                                 ignore.case  = FALSE, 
                                 include.dirs = FALSE) %>%
                    Rsamtools::BamFileList(
                                         yieldSize = 7500000)

register( MulticoreParam( workers = 6) )
ATACseq.counts           <- summarizeOverlaps(  
                                  features      = common_Granges,
                                  reads         = ATAC.peak.bams, 
                                  ignore.strand = TRUE, 
                                  singleEnd     = TRUE ) %>% assay()

# new start
#---

ATACseq.mm10.annot     <- toGRanges( TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                   feature = 'gene')

seqlevelsStyle(common_Granges) <- seqlevelsStyle(ATACseq.mm10.annot)
ATACseq.normed.counts <- ATACseq.counts %>%
                        cpm(normalized.lib.sizes = TRUE, log = TRUE)

ATACseq.peaks.annot        <- common_Granges %>%
                            annotatePeakInBatch( ,
                                   AnnotationData = ATACseq.mm10.annot,
                                   select         = 'first')

ATACseq.symbols  <- mapIds( org.Mm.eg.db, 
                                    keys      = ATACseq.peaks.annot  %>%
                                                mcols() %$% feature,
                                    column    = 'SYMBOL', 
                                    keytype   = 'ENTREZID', 
                                    multiVals = 'first') %>%
                            make.names(unique = TRUE)
rownames(ATACseq.counts) <- ATACseq.symbols
colnames(ATACseq.counts) <- c( 'P1_1', 'P1_2', 'P1_3',
                               'P14_1','P14_2','P14_3',
                               'P56_1','P56_2','P56_3') 

color.bar     <- colorRampPalette(c('blue4', 'white', 'springgreen4'))(10) 
color.bar     <- colorRampPalette(brewer.pal(10,'RdYlBu'))(10) %>% rev()
cib1.pheatmap <- ATACseq.counts  %>% {cor(., method = 'pearson')} %>%
                 pheatmap( cluster_rows = T, cluster_cols = F,
                           color = c(color.bar), scale = 'none', fontsize_col = 8, 
                           fontsize_row = 8, cellwidth = 30, cellheight = 30,
                           labels_row = colnames(ATACseq.counts), 
                           labels_col = colnames(ATACseq.counts),
                           display_numbers = TRUE, number_color = 'orange',
                           fontsize_number = 10)

cib1.pheatmap <- ATACseq.counts   %>%
                 pheatmap( cluster_rows = T, cluster_cols = F,
                           color = c(color.bar), scale = 'row',  
                           labels_col = colnames(ATACseq.counts),
                           labels_row = FALSE, show_rownames = FALSE,
                           cutree_rows = 5)

fviz_nbclust(ATACseq.counts, pam, method = "silhouette")+
theme_classic()

hclust(dist(ATACseq.counts))

pheatmap(ATACseq.counts, cutree_row = 3)

save.image('ATACseq.Rdata')
