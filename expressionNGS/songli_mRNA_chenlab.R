# @author Yisong Zhen
# @since  2018-04-20
# @update 2018-04-23
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


mRNA.songli.path    <- file.path('/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs')
hg38.rusbread.index <- file.path('/home/zhenyisong/data/reference/index')
output.path         <- file.path('/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs/rsubread')
QC.path             <- file.path('/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs/rsubread/multiQC')

ifelse(!dir.exists(output.path), dir.create(output.path), FALSE)
ifelse(!dir.exists(QC.path), dir.create(QC.path), FALSE)

"
working.env <- 'window'
linux.path  <- mRNA.songli.path 
window.path <- file.path('D:\\yisong.data')
image.data  <- 'mRNA.songli.Rdata'


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
cluster         <- makeCluster(3)

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
save.image('mRNA.songli.QC.Rdata')
q('no')

setwd(hg38.rusbread.index)

base.string        <-  'hg38'  
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

mRNA.count.table <- mRNA.songli.genes %$% counts
colnames(mRNA.count.table) <- c( 'SY001','SY037','SY045',
                                 'SY052','SY055','SY057',
                                 'SY062','SY065','SY086',
                                 'SY091','SY097','SY111')   
mRNA.rsubread.QC           <- mRNA.songli.genes$stat
colnames(mRNA.rsubread.QC) <- c('status', colnames(mRNA.count.table))

"
In songli's email coorespondings   Group B: high expression hi_exprs                                    
Group A: low exression lo_exprs    SY111                               
SY086                              SY055     
SY062                              SY045     
SY052                              SY001     
SY065                              SY037     
SY091                              SY097     
                                   SY057
"



mRNA.group       <- factor( c( 1, 1, 1, 0, 1, 1,
                               0, 0, 0, 0, 1, 1  ), 
                            levels = 0:1,
                            labels = c( 'lo_exprs','hi_exprs'))
column.mRNA           <- data.frame(design = mRNA.group)
rownames(column.mRNA) <- colnames(mRNA.count.table)
diff.deseq2    <- DESeqDataSetFromMatrix( countData = mRNA.count.table,
                                                 colData   = column.mRNA,
                                                 design    = ~ design ) %>%
                  DESeq()
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
                     mutate(symobls = gene.names) %>%
                     arrange(pvalue)

songli.mRNA.wb <- createWorkbook()
addWorksheet(songli.mRNA.wb, 'mRNA_diff')
addWorksheet(songli.mRNA.wb, 'mRNA_rusubread_QC')
writeData(songli.mRNA.wb, sheet = 1, diff.output.table )
writeData(songli.mRNA.wb, sheet = 2, mRNA.rsubread.QC )
saveWorkbook(songli.mRNA.wb, 'songlimRNA.2018-04-24.xlsx', overwrite = TRUE)


setwd(output.path)
save.image('mRNA.songli.Rdata')