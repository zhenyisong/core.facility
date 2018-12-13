# Title     : Her Huamn patient sample QC assay
# Objective : Genenrate PCA post-processing plot to be publihsed
# Created by: zhenyisong
# Created on: 12/13/2018
#
#----

#---
# now generate the fingerprints of all the RNAseq data for validaition
#---

"
find -type f -name '*.fq.gz' -exec md5sum '{}' + > chenlab_data.fingerprints
"


pkgs     <- c( 'tidyverse','Rsubread','org.Hs.eg.db','edgeR',
               'limma', 'DESeq2', 'genefilter','grid',
               'openxlsx','pheatmap','gridExtra','ggrepel',
               'QuasR','annotate','clusterProfiler',
               'cowplot','magrittr', 'reshape2',
               'GGally','RColorBrewer', 'PROPER',
               'cluster','factoextra','ggpubr',
               'BSgenome.Hsapiens.UCSC.hg38',
               'BSgenome.Hsapiens.UCSC.hg38.Rbowtie')

load.lib            <- lapply(pkgs, require, character.only = TRUE)

hg38.rusbread.index <- file.path('/wa/zhenyisong/reference/index')
raw.data.path       <- file.path('/wa/zhenyisong/results/chenlab/upload/Cleandata')
output.bam.path     <- file.path('/wa/zhenyisong/results/chenlab/upload/rsubread')
output.qc.path      <- file.path('/wa/zhenyisong/results/chenlab/upload/rsubread/multiQC')

genome.specie       <-  'BSgenome.Hsapiens.UCSC.hg38'

"
ifelse(!dir.exists(output.bam.path), dir.create(output.bam.path), FALSE)
ifelse(!dir.exists(output.qc.path), dir.create(output.qc.path), FALSE)
"

image.data  <- 'qc.chenxi.Rdata'


get.complete.rawdata.list <- function(data.path, suffix.pattern) {
    raw.data <- list.files( path         = data.path,
                            pattern      = suffix.pattern,
                            all.files    = FALSE,
                            full.names   = TRUE,
                            recursive    = TRUE,
                            ignore.case  = FALSE,
                            include.dirs = TRUE)
    return(raw.data)

}

seperate.PE.rawdata   <- function(rawdata.list, pattern.1 = 'R1', pattern.2 = 'R2') {
    read.1.files   <- grep(pattern.1,rawdata.list) %>% rawdata.list[.]
    read.2.files   <- grep(pattern.2,rawdata.list) %>% rawdata.list[.]
    seperated.list <- list('read.1' = read.1.files, 'read.2' = read.2.files)
    return(seperated.list)
}

get.sample.names <- function(read.1.files, pattern = '_R1.fq.gz') {
    sample.names    <- basename(read.1.files) %>%
                       sub(pattern = pattern, replacement = '')
    return(sample.names)
}


get.ouptut.file.names <- function(sample.names,output.path) {
    output.file.names <- paste(output.path,'/', sample.names,'.bam', sep = '')
    return(output.file.names)
}

.get.temp.sampleFile <- function(read.1.files, read.2.files, sampleNames) {
    sampleFile              <- tempfile(pattern = 'zhen3temp', tmpdir = tempdir())
    sample.file             <- data.frame( FileName1  = read.1.files,
                                           FileName2  = read.2.files,
                                           SampleName = sampleNames)
    write_tsv(sample.file, path = sampleFile)
    return(sampleFile)
}


.set_QC_check_parameters <- function(specie = genome.specie,
                                     cpus  = 5, qc.result.path) {
    setwd(qc.result.path)
    genome     <- specie
    cluster    <- makeCluster(cpus)
    param.list <- list('genome' = genome, 'cluster' = cluster, 'output' = qc.result.path)
    return(param.list)
}

get_qc_final_obj <- function( file_pattern = 'fq.gz$',
                              read.1 = 'R1', read.2 = 'R2',
                              read.1.suffix = '_R1.fq.gz',
                              specie = genome.specie,
                              cpus   = 5,
                              full.data.path  =  raw.data.path,
                              qc.output.path = output.qc.path) {
    qc.params    <- .set_QC_check_parameters(qc.result.path = qc.output.path)
    rawdata.list <- get.complete.rawdata.list(full.data.path, file_pattern)
    reads.list   <- seperate.PE.rawdata(rawdata.list,read.1, read.2)
    sample.names <- get.sample.names(reads.list$read.1, read.1.suffix)
    sampleFile   <- .get.temp.sampleFile(reads.list$read.1,reads.list$read.2,sample.names)
    #print(qc.params)
    #print(rawdata.list)
    #print(sampleFile)
    #print(reads.list)
    QC.qPorject  <- qAlign( sampleFile,
                            qc.params$genome,
                            auxiliaryFile = NULL,
                            aligner = 'Rbowtie',
                            maxHits = 1,
                            paired  = 'fr',
                            splicedAlignment = FALSE,
                            snpFile   = NULL,
                            bisulfite = 'no',
                            alignmentParameter = NULL,
                            projectName   = 'corefacility',
                            alignmentsDir = qc.params$output,
                            lib.loc   = NULL,
                            cacheDir  = NULL,
                            clObj     = qc.params$cluster,
                            checkOnly = F)
    qQCReport( QC.qPorject,
               pdfFilename    = 'final.QC.pdf',
               useSampleNames = TRUE,
               clObj          = qc.params$cluster)
    stopCluster(qc.params$cluster)
}


get_rsubread_result <- function(genome.index.path = hg38.rusbread.index ,
                                index.tag         = 'hg38',
                                rsubread.genome   = 'hg38',
                                full.data.path    =  raw.data.path,
                                file_pattern = 'fq.gz$',
                                read.1 = 'R1', read.2 = 'R2',
                                read.1.suffix = '_R1.fq.gz',
                                output.path   = output.bam.path
                                ) {
    setwd(genome.index.path)
    rawdata.list        <- get.complete.rawdata.list(full.data.path, file_pattern)
    reads.list          <- seperate.PE.rawdata(rawdata.list,read.1, read.2)
    sample.names        <- get.sample.names(reads.list$read.1, read.1.suffix)
    output.filenames    <- get.ouptut.file.names(sample.names,output.path)
    base.string         <- index.tag
    #print(output.filenames)

    align.result        <- align( index          = base.string,
                                  readfile1      = reads.list$read.1,
                                  readfile2      = reads.list$read.2,
                                  input_format   = 'gzFASTQ',
                                  type           = 'rna',
                                  output_file    = output.filenames,
                                  output_format  = 'BAM',
                                  PE_orientation = 'fr',
                                  nthreads       = 6,
                                  indels         = 1,
                                  maxMismatches  = 3,
                                  phredOffset    = 33,
                                  unique         = T )
    rsubread.result     <- featureCounts( output.filenames,
                                          useMetaFeatures        = TRUE,
                                          countMultiMappingReads = FALSE,
                                          strandSpecific         = 0,
                                          isPairedEnd            = TRUE,
                                          requireBothEndsMapped  = TRUE,
                                          autosort               = TRUE,
                                          nthreads               = 6,
                                          annot.inbuilt          = rsubread.genome,
                                          allowMultiOverlap      = TRUE)
    return(rsubread.result)
}


get_power_analysis_result <- function(gene.counts) {

    phenotype.param  <- estParam(gene.counts)

    "
    The following command sets up simulation options with 20,000 genes, 5% genes being DE,
    baseline expression and dispersion based on Chenlab's data
    "

    sim.opt.human    <- RNAseq.SimOptions.2grp( ngenes = 20000,
                                                p.DE   = 0.05,
                                                lOD    = phenotype.param$lOD,
                                                lBaselineExpr = phenotype.param$lmean)
    "
    By default, we evaluate power when there are 16 ~ 112
    samples in each treatment group
    "
    phenotype.simres <- runSims( Nreps    = c(16,32,48,64),
                                 sim.opts = sim.opt.human,
                                 DEmethod = 'edgeR',
                                 nsims    = 100)
    return(phenotype.simres)
}

#---
# https://github.com/zhenyisong/wanglab.code/blob/master/yaoyan.R
#---

set_rownames_as_genesymbol <- function(vsd.exprs, rsubread.result) {
    gene.ids            <- rsubread.result$annotation$GeneID
    gene.info           <- AnnotationDbi::select(
                               org.Hs.eg.db, keys= as.character(gene.ids),
                               keytype = 'ENTREZID', columns = 'SYMBOL')
    rownames(vsd.exprs) <- gene.info$SYMBOL %>% make.names(unique = T)
    return(vsd.exprs)
}


get_sex_determination_pca_result <- function(vsd.exprs) {
    mart         <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
    results      <- getBM( attributes = c('chromosome_name', 'entrezgene', 'hgnc_symbol'),
                           filters    = 'chromosome_name', values = 'Y', mart = mart)
    y.chromosome.gene <- unique(results$hgnc_symbol) %>% na.omit()
    y.chromosome.gene <- y.chromosome.gene[ y.chromosome.gene != '']
    vsd.expr.sex      <- vsd.exprs[ rownames(vsd.exprs) %in% y.chromosome.gene,]
    #vsd.expr.sex      <- vsd.pca[ rownames(vsd.pca) %in% y.chromosome.gene,]
    #colnames(vsd.expr.sex)
    #colnames(vsd.expr.sex) <- sample.names
    sds.sex <- rowSds(vsd.expr.sex)
    summary(sds.sex)
    vsd.expr.sex[sds.sex > 0.5,c(1,6,7,11)]
    km.res <- kmeans(t(vsd.expr.sex[sds.sex > .7,]), 2, nstart = 100)


    nu.number <- fviz_nbclust( t(vsd.expr.sex[sds.sex > .7,]),
                               kmeans, method = "silhouette") +
                 theme_classic()

    cluster.figure <- fviz_cluster(km.res, data = t(vsd.expr.sex[sds.sex > .7,]),
                                   palette = c('#2E9FDF', '#FC4E07'),
                                   ellipse.type = 'euclid', # Concentration ellipse
                                   star.plot    = TRUE, # Add segments from centroids to items
                                   repel        = TRUE, # Avoid label overplotting (slow)
                                   ggtheme      = theme_minimal())
    result <- list('cluster.number' = nu.number, 'full.figure' = cluster.figure)
    return(result)

}


get_rsubread_qc_result <- function(rsubread.result,colunm.names) {
    rsubread.qc   <- rsubread.result$stat
    colnames(rsubread.qc) <- c('status', colunm.names)
    return(rsubread.qc)
}


get_vsd_for_pca_analysis <- function(rsubread.results,column.names, groups.design) {
    mRNA.count.table           <- rsubread.results %$% counts
    colnames(mRNA.count.table) <- complete.sample.names
    column.mRNA                <- data.frame(design = groups.design)
    rownames(column.mRNA)      <- column.names
    diff.deseq2                <- DESeqDataSetFromMatrix( countData = mRNA.count.table,
                                                          colData   = column.mRNA,
                                                          design    = ~ design) %>%
                                  DESeq()
    mRNA.vsd                   <- getVarianceStabilizedData(diff.deseq2)
    return(mRNA.vsd)
}


get_pca_cluster <- function(vsd.result,groups, sample.names) {
    mRNA.pca           <- prcomp(t(vsd.result))
    pca.no             <- dim(mRNA.pca$x)[2]
    graph.text.size    <- 8
    point.size         <- 2
    mRNA.PCA.ggplot    <- as.data.frame(mRNA.pca$x) %>%
                          mutate(color.choice = groups) %>%
                          ggplot(data = . ) +
                          geom_point(aes(x = PC1, y = PC2, color = color.choice), size = point.size) +
                          geom_text_repel(aes(x = PC1, y = PC2, label = sample.names)) +
                          scale_colour_manual( name   = 'group in experiment design',
                                               values = c( 'darkblue','green','red'),
                                               labels = c( 'control','border','patient')) +
                          theme_bw() +
                          theme(legend.position = c(.5,.5)) +
                          theme( aspect.ratio = 1, text = element_text(size = graph.text.size),
                                 legend.text  = element_text(size = 12),
                                 legend.title = element_text(size = 12))
    return(mRNA.PCA.ggplot)
}


get_pca_cumsum_plot <- function(vsd.result,groups, sample.names) {
    mRNA.pca           <- prcomp(t(vsd.result))
    pca.no             <- dim(mRNA.pca$x)[2]
    graph.text.size    <- 8
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
                          theme_bw() +
                          theme( aspect.ratio = 1, text = element_text(size = graph.text.size))
    return(cumsum.plot.ggplot)
}


#setwd(output.bam.path)
#ene.counts  <- get_rsubread_result()

#setwd(output.bam.path)
setwd('C:\\Users\\zheny\\Downloads')
load(image.data)

"
the sample coding from Chen's student who grouped
three kinds of samples together
"

control_samples <- c('11-L', '223P', '371P',
                     '392L', '483', '588', '594', '605')
control_samples <- paste0('control_', control_samples)
border_samples  <- c('106', '131', '34-L',
                     '404P', '455', '525', '56L','596')
border_samples  <- paste0('border_', border_samples)
patient_samples <- c('126', '142', '143',
                     '160', '208P', '37-L', '38L','461')
patient_samples <- paste0('patient_',patient_samples)


human.sample.design      <- factor( c(rep(0,8),rep(1,8),rep(2,8)),
                                    levels = 0:2,
                                    labels = c('control','border','patient'))
complete.sample.names    <- c(control_samples, border_samples,patient_samples)
rsubread.qc              <- get_rsubread_qc_result(rsubread.result,complete.sample.names)

vsd.pca                  <- get_vsd_for_pca_analysis( rsubread.result,
                                                      complete.sample.names,
                                                      human.sample.design)
vsd.pca          <- set_rownames_as_genesymbol(vsd.pca, rsubread.result)
pr.out           <- prcomp(t(vsd.pca))
pr.var           <- pr.out$sdev^2
pve              <- pr.var/sum(pr.var)

sex.3.figure     <- get_sex_determination_pca_result(vsd.pca)
pca.1.figure     <- get_pca_cumsum_plot(vsd.pca, human.sample.design,complete.sample.names)
pca.2.figure     <- get_pca_cluster(vsd.pca, human.sample.design,complete.sample.names)
power.simulation <- get_power_analysis_result(rsubread.result$counts[,c(1:8,17:24)])

phenotype.power  <- comparePower(power.simulation)
summaryPower(phenotype.power)
phenotype.power.update <- comparePower( power.simulation,
                                        alpha.nominal = 0.00001,
                                        alpha.type    = 'pval',
                                        delta         = log(1))
power.update <- summaryPower(phenotype.power.update)
plotPower(phenotype.power.update)
plotAll(phenotype.power.update)
g <- tableGrob(power.update, rows = NULL)
grid.newpage()
grid.draw(g)

"
save.image(image.data)
get_qc_final_obj()
"
save.image('power_analysis.Rdata')







