# @author Yisong Zhen
# @since  2018-05-04
# @update 2018-05-08
#---

pkgs <- c( 'tidyverse','Rsubread','org.Hs.eg.db','edgeR',
           'limma', 'DESeq2', 'genefilter','grid',
           'openxlsx','pheatmap','gridExtra','ggrepel',
           'QuasR','annotate','clusterProfiler',
           'cowplot','magrittr',
           'GGally','RColorBrewer',
           'cluster','factoextra','ggpubr',
           'BSgenome.Hsapiens.UCSC.hg38',
           'BSgenome.Hsapiens.UCSC.hg38.Rbowtie')

load.lib <- lapply(pkgs, require, character.only = TRUE)


mRNA.songli.path    <- file.path('/wa/zhenyisong/results/chenlab/songli/phenotype')
hg38.rusbread.index <- file.path('/home/zhenyisong/data/reference/index')
output.path         <- file.path('/wa/zhenyisong/results/chenlab/songli/phenotype/rsubread')
QC.path             <- file.path('/wa/zhenyisong/results/chenlab/songli/phenotype/rsubread/multiQC')

ifelse(!dir.exists(output.path), dir.create(output.path), FALSE)
ifelse(!dir.exists(QC.path), dir.create(QC.path), FALSE)

"
working.env <- 'window'
linux.path  <- mRNA.songli.path 
window.path <- file.path('D:\\yisong.data')
image.data  <- 'phenotype.songli.Rdata'


# rm(list=ls())

switch( working.env, linux  = { setwd(linux.path);
                                load( file = image.data)},
                     window = { setwd(window.path);
                                load( file = image.data)} )
"

mRNA.songli.data        <- list.files( path         = mRNA.songli.path, 
                                       pattern      = '.clean.fq.gz$', 
                                       all.files    = FALSE, 
                                       full.names   = TRUE, 
                                       recursive    = FALSE, 
                                       ignore.case  = FALSE, 
                                       include.dirs = FALSE)
   
read.1.files            <- grep('_1',mRNA.songli.data) %>% 
                           mRNA.songli.data[.]
read.2.files            <- grep('_2',mRNA.songli.data) %>% 
                           mRNA.songli.data[.]

mRNA.sample.names       <- basename(read.1.files) %>% 
                           sub(pattern = '_1.clean.fq.gz', replacement = '')
mRNA.output.filenames   <- paste(output.path,'/', mRNA.sample.names,'.bam', sep = '')

sampleFile              <- tempfile(pattern = 'zhen3temp', tmpdir = tempdir())
sample.file             <- data.frame( FileName1  = read.1.files,
                                       FileName2  = read.2.files, 
                                       SampleName = mRNA.sample.names)

write_tsv(sample.file, path = sampleFile)
genome          <- 'BSgenome.Hsapiens.UCSC.hg38'
cluster         <- makeCluster(5)

setwd(QC.path)

songli.qPorject  <- qAlign( sampleFile,
                            genome,
                            auxiliaryFile = NULL,
                            aligner = 'Rbowtie',
                            maxHits = 1,
                            paired  = 'fr',
                            splicedAlignment = FALSE,
                            snpFile = NULL,
                            bisulfite = 'no',
                            alignmentParameter = NULL,
                            projectName = 'songli',
                            alignmentsDir = QC.path,
                            lib.loc  = NULL,
                            cacheDir = NULL,
                            clObj = cluster,
                            checkOnly = F)

qQCReport( songli.qPorject, 
           pdfFilename    = 'songli.QC.pdf', 
           useSampleNames = TRUE, 
           clObj          = cluster)
stopCluster(cluster)

#save.image('mRNA.songli.QC.Rdata')
#q('no')

setwd(hg38.rusbread.index)

base.string        <- 'hg38'  
align.result       <- align( index          = base.string, 
                             readfile1      = read.1.files, 
                             readfile2      = read.2.files, 
                             input_format   = 'gzFASTQ', 
                             type           = 'rna',
                             output_file    = mRNA.output.filenames, 
                             output_format  = 'BAM',
                             PE_orientation = 'fr', 
                             nthreads       = 6, 
                             indels         = 1,
                             maxMismatches  = 3,
                             phredOffset    = 33,
                             unique         = T )
mRNA.songli.genes  <- featureCounts( mRNA.output.filenames, 
                             useMetaFeatures        = TRUE,
                             countMultiMappingReads = FALSE,
                             strandSpecific         = 0, 
                             isPairedEnd            = TRUE,
                             requireBothEndsMapped  = TRUE,
                             autosort               = TRUE,
                             nthreads               = 6,
                             annot.inbuilt          = 'hg38', 
                             allowMultiOverlap      = TRUE)

"
setwd(output.path)
save.image('phenotype.songli.Rdata')
q('no')
"
get_sample_names <- function( rsubread.aligned, pattern = '\\.(\\w+)\\.bam') {
    sample.names   <- rsubread.aligned %$% 
                      counts %>%
                      colnames() %>%
                      map(str_match, pattern) %>%
                      map(2) %>%
                      unlist()
    return(sample.names)
}

get_Rsubread_QC  <- function(rsubread.aligned, pattern = '\\.(\\w+)\\.bam') {
    rsubread.QC           <- rsubread.aligned$stat
    sample.names          <- get_sample_names(rsubread.aligned, pattern)
    colnames(rsubread.QC) <- c('status', sample.names)
    return(rsubread.QC)
}

phenotype.QC.table <- get_Rsubread_QC(mRNA.songli.genes)
(phenotype.QC.table)

groups <- factor( c( 0, 0, 1, 1, 0, 1,
                     1, 1, 0, 0, 1, 0,
                     0, 1, 0, 1, 0, 1,
                     1 ), 
                  levels = 0:1,
                  labels = c( 'snp_GG','snp_TT'))

get_DESeq2_object <- function( groups.factor, rsubread.aligned, pattern) {
    column.names            <- data.frame(design = groups.factor)
    rownames(column.names)  <- get_sample_names(rsubread.aligned, pattern)
    counts.matrix           <- rsubread.aligned %$% counts
    colnames(counts.matrix) <- rownames(column.names)
    diff.deseq2             <- DESeqDataSetFromMatrix( 
                                    countData = counts.matrix,
                                    colData   = column.names,
                                    design    = ~ design ) %>%
                               DESeq()
    return(diff.deseq2)        
}

get_DESeq2_object(groups, mRNA.songli.genes, pattern = '\\.(\\w+)\\.bam') %>%
getVarianceStabilizedData()


mRNA.vsd       <- getVarianceStabilizedData(diff.deseq2)
heatmap(cor(mRNA.vsd))
mRNA.pca <- prcomp(t(mRNA.vsd)) 
plot( mRNA.pca$x, 
      col  = 'white', 
      main = 'PCA plot', 
      xlim = c(-70,150),
      ylim = c(-40,110))
text( mRNA.pca$x[,1], mRNA.pca$x[,2], 
      labels = colnames(mRNA.vsd),
      cex    = 0.6)

pca.no             <- dim(mRNA.pca$x)[2]
graph.text.size    <- 4
point.size         <- 2
cumsum.plot.ggplot <- {mRNA.pca$sdev^2} %>%
                      {./sum(.)} %>% 
                      {data.frame(variance = cumsum(.), pca = c(1:pca.no))}%>%
                      ggplot(data = .) +
                      xlab('Principle Components') +
                      ylab('Culmulative Variance') +
                      scale_x_continuous( breaks = c(1:pca.no), labels = as.character(c(1:pca.no), 
                                          limits = as.character(c(1:pca.no)))) +
                      geom_point(aes(x = pca, y = variance), size = point.size) +
                      geom_line(aes(x = pca, y = variance),size = 0.8) +
                      scale_linetype_discrete() +
                      theme(legend.position = 'none') +
                      theme_classic() +
                      theme( aspect.ratio = 1, text = element_text(size = graph.text.size))
(cumsum.plot.ggplot)


mRNA.PCA.ggplot <- as.data.frame(mRNA.pca$x) %>% 
                   mutate(color.choice = mRNA.group) %>%
                   ggplot(data = . ) + 
                   geom_point(aes(x = PC1, y = PC2, color = color.choice), size = point.size) + 
                   geom_text_repel(aes(x = PC1, y = PC2, label = rownames(column.mRNA))) +
                   scale_colour_manual( name   = 'group in experiment design',
                                        values = c( 'darkblue','green'),
                                        labels = c( 'low_exprs','high_exprs')) +   
                   theme_classic() +
                   theme( aspect.ratio = 1, text = element_text(size = graph.text.size))
(mRNA.PCA.ggplot)

comp.name   <- resultsNames(diff.deseq2)[2]
diff.result <- results(diff.deseq2, name = comp.name)

gene.names  <- mapIds( org.Hs.eg.db, 
                       keys      = rownames(diff.result) %>% 
                                   as.character(),
                       column    = 'SYMBOL', 
                       keytype   = 'ENTREZID', 
                       multiVals = 'first') %>%
                       make.names(unique = T)


diff.output.table <- diff.result %>% 
                     as.data.frame %>% 
                     mutate(symbols = gene.names) %>%
                     arrange(pvalue)

# extra steps to combine the ChIRPseq and mRNAseq results
# we solicit the common genes in both sets.
# 'peaks.hg38.gene.symbol' object is from 
# blood_ChIRP_songli.R
# and import from 'songliChIRP.Rdata'.
#---
mRNA.gene.names <- diff.output.table %>% 
                   filter(pvalue <= 0.05) %>%
                   dplyr::select(symbols) %>% unlist()

ChIRP.gene.names <- peaks.hg38.gene.symbols %>% unlist()

common.gene.names <- intersect(mRNA.gene.names, ChIRP.gene.names)

songli.mRNA.wb <- createWorkbook()
addWorksheet(songli.mRNA.wb, 'mRNA_diff')
addWorksheet(songli.mRNA.wb, 'mRNA_rusubread_QC')
addWorksheet(songli.mRNA.wb, 'two_set_common_genes')
writeData(songli.mRNA.wb, sheet = 1, diff.output.table )
writeData(songli.mRNA.wb, sheet = 2, mRNA.rsubread.QC )
writeData(songli.mRNA.wb, sheet = 3, common.gene.names )
saveWorkbook(songli.mRNA.wb, 'songlimRNA.2018-04-24.xlsx', overwrite = TRUE)



setwd(output.path)
save.image('mRNA.songli.Rdata')