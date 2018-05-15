# @author Yisong Zhen
# @since  2018-05-04
# @update 2018-05-11
#---

pkgs     <- c( 'tidyverse','Rsubread','org.Hs.eg.db','edgeR',
               'limma', 'DESeq2', 'genefilter','grid',
               'openxlsx','pheatmap','gridExtra','ggrepel',
               'QuasR','annotate','clusterProfiler',
               'cowplot','magrittr', 'reshape2',
               'GGally','RColorBrewer', 'PROPER',
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


working.env <- 'window'
linux.path  <- output.path 
window.path <- file.path('D:\\yisong.data')
image.data  <- 'phenotype.songli.Rdata'


# rm(list=ls())

switch( working.env, linux  = { setwd(linux.path);
                                load( file = image.data)},
                     window = { setwd(window.path);
                                load( file = image.data)} )


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


#---
# power calculation
#---
phenotype.param  <- estParam(mRNA.songli.genes$counts)
sim.opt.human    <- RNAseq.SimOptions.2grp( ngenes = 20000, 
                                            p.DE   = 0.05, 
                                            lOD    = phenotype.param$lOD,
                                            lBaselineExpr = phenotype.param$lmean)
#simRNAseq(simOptions, n1 = 3, n2 = 3)
   
phenotype.simres <- runSims( Nreps    = c(9, 12, 15), 
                             sim.opts = sim.opt.human, 
                             DEmethod = 'edgeR', 
                             nsims    = 100)
phenotype.power  <- comparePower(phenotype.simres)
summaryPower(phenotype.power)
phenotype.power.update <- comparePower( phenotype.simres, 
                                        alpha.nominal = 0.05, 
                                        alpha.type    = 'pval',
                                        delta         = log(0.5))
power.simu.results <- summaryPower(phenotype.power.update)
grid.newpage() %>% 
{tableGrob( power.simu.results, 
            rows = NULL)}

g <- tableGrob(power.simu.results, rows = NULL)
grid.newpage()
grid.draw(g)

plotAll(phenotype.power.update)




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



get_transformed_vsd <- function(catergories, rusbread.aligned, pattern) {
    data.vsd <- get_DESeq2_object( catergories, 
                                   rusbread.aligned, 
                                   pattern = pattern) %>%
                getVarianceStabilizedData()
    return(data.vsd)
}


get_pca_data <- function(vsd_trasnformed_data) {
    data.pca <- prcomp(t(vsd_trasnformed_data)) 
    return(data.pca)
}

reorder_cormat <- function(data) {
    hc     <- as.dist((1 - data)/2) %>%
              hclust(.)
    cormat <- data[hc$order, hc$order]
    return(data)
}

get_lower_tri <- function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}

get_postQC_figures <- function( rsubread.aligned,
                                pattern, 
                                catergories, 
                                group.colors,
                                group.labels) {
    sample.names <- get_sample_names(rsubread.aligned, pattern)
    data.vsd   <- get_transformed_vsd( catergories, 
                                       rsubread.aligned,
                                       pattern = pattern )
    pca.data   <- get_pca_data(data.vsd)
    pca.no             <- dim(pca.data$x)[2]
    graph.text.size    <- 4
    point.size         <- 2
    cumsum.plot.ggplot <- {pca.data$sdev^2} %>%
                          {./sum(.)} %>% 
                          {data.frame( variance = cumsum(.), 
                                       pca = c(1:pca.no)) } %>%
                          ggplot(data = .) +
                          xlab('Principle Components') +
                          ylab('Culmulative Variance') +
                          scale_x_continuous( breaks = c(1:pca.no), 
                                              labels = as.character(c(1:pca.no), 
                                              limits = as.character(c(1:pca.no)))) +
                          geom_point( aes( x   = pca, 
                                           y   = variance), 
                                      size = point.size ) +
                          geom_line( aes( x    = pca, 
                                          y    = variance),
                                     size = 0.8) +
                          scale_linetype_discrete() +
                          theme(legend.position = 'none') +
                          theme_classic() +
                          theme( aspect.ratio = 1, 
                                 text = element_text( size = graph.text.size))
    mRNA.PCA.ggplot <- as.data.frame(pca.data$x) %>% 
                       mutate(color.choice = catergories) %>%
                       ggplot(data = . ) + 
                       geom_point( aes( x = PC1, 
                                        y = PC2, 
                                        color = color.choice), 
                                   size = point.size) + 
                       geom_text_repel( aes( x = PC1, 
                                             y = PC2, 
                                             label = sample.names)) +
                       scale_colour_manual( name   = 'group in experiment design',
                                            values = group.colors,
                                            labels = group.labels)  +
                       theme_classic() +
                       theme( aspect.ratio = 1, legend.position = 'none',
                              text = element_text(size = graph.text.size))

    ggheatmap    <- round(cor(data.vsd),2) %>%
                     reorder_cormat() %>%
                     get_upper_tri() %>%
                     melt(na.rm = TRUE) %>%
                     ggplot(aes(Var2, Var1, fill = value)) +
                     geom_tile(color = "white") +
                     scale_fill_gradient2( low  = "blue", 
                                           high = "red", 
                                           mid = "white", 
                                           midpoint = 0, 
                                           limit = c(-1,1), 
                                           space = "Lab", 
                                           name  = "Pearson\nCorrelation") +
                     theme_minimal() + # minimal theme
                     theme( axis.text.x = 
                                      element_text( angle = 45, vjust = 1, 
                                                    size = 12, hjust = 1),
                            aspect.ratio = 1) +
                     coord_fixed()

    scree.ggplot <-  {pca.data$sdev^2} %>%
                     {./sum(.)} %>% 
                     {data.frame(variance = ., pca = c(1:pca.no)) } %>%
                     ggplot(data = .) +
                     xlab('Wangyin lincRNA Principle Components') +
                     ylab('Proportion of Variance Explained') +
                     scale_x_continuous( breaks = c(1:pca.no), 
                                         labels = as.character(c(1:pca.no), 
                                         limits = as.character(c(1:pca.no)))) +
                     geom_point(aes(x = pca, y = variance), size = 3) +
                     geom_line(aes(x = pca, y = variance), size = 0.8) +
                     scale_linetype_discrete() +
                     theme_classic() +
                     theme(legend.position = 'none', aspect.ratio = 1) 
                     
    final.figures <- ggdraw() +
                     draw_plot(ggheatmap,       x =  0, y = .5,  width = .5, height = .5) +
                     draw_plot(scree.ggplot,    x = .5, y = .5,  width = .5, height = .5) +
                     draw_plot(cumsum.plot.ggplot,   x =  0, y =  0,  width = .5, height = .5) +
                     draw_plot(mRNA.PCA.ggplot,  x = .5, y =  0,  width = .5, height = .5) +
                     draw_plot_label( label = LETTERS[1:4],
                                      x     = c(0,.5, 0,.5),
                                      y     = c(1, 1,.5,.5),
                                      size  = 10)
    return(final.figures)
}


get_diff_genes <- function(DGE.deseq2, org = org.Hs.eg.db) {

    gene.names        <- mapIds( org, 
                                 keys      = rownames(DGE.deseq2) %>% 
                                             as.character(),
                                 column    = 'SYMBOL', 
                                 keytype   = 'ENTREZID', 
                                 multiVals = 'first') %>%
                                 make.names(unique = T)
      

    diff.output.table <- DGE.deseq2 %>% 
                         as.data.frame %>% 
                         mutate(symbols = gene.names) %>%
                         arrange(pvalue)
    return(diff.output.table)
}




phenotype.QC.table <- get_Rsubread_QC(mRNA.songli.genes)
(phenotype.QC.table)

groups           <- factor( c( 0, 0, 1, 1, 0, 1,
                               1, 1, 0, 0, 1, 0,
                               0, 1, 0, 1, 0, 1,
                               1 ), 
                            levels = 0:1,
                            labels = c( 'snp_GG','snp_TT'))
          
postQC.figures   <- get_postQC_figures( mRNA.songli.genes,
                                        pattern = '\\.(\\w+)\\.bam', 
                                        groups, 
                                        group.colors = c('red','blue'),
                                        group.labels = c('GG','TT'))

diff.deseq2      <- get_DESeq2_object( groups, mRNA.songli.genes, 
                                       pattern = '\\.(\\w+)\\.bam')
comp.name        <- resultsNames(diff.deseq2)[2]
diff.result      <- results(diff.deseq2, name = comp.name)

final.DGE.result <- get_diff_genes(diff.result, org.Hs.eg.db)


songli.phenotype.wb <- createWorkbook()
addWorksheet(songli.phenotype.wb, 'mRNA_diff')
addWorksheet(songli.phenotype.wb, 'mRNA_rusubread_QC')
writeData(songli.phenotype.wb, sheet = 1, final.DGE.result )
writeData(songli.phenotype.wb, sheet = 2, phenotype.QC.table )
saveWorkbook(songli.phenotype.wb, 'songliPhenotype.2018-05-11.xlsx', overwrite = TRUE)



setwd(output.path)
save.image('mRNA.songli.Rdata')