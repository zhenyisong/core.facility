# Title     : Her Huamn patient sample QC assay
# Objective : Genenrate PCA post-processing plot to be publihsed
# Created by: zhenyisong
# Created on: 12/07/2018
#----


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

ifelse(!dir.exists(output.bam.path), dir.create(output.bam.path), FALSE)
ifelse(!dir.exists(output.qc.path, dir.create(output.qc.path), FALSE)

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
    print(qc.params)
    print(rawdata.list)
    print(sampleFile)
    print(reads.list)
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
    print(output.filenames)

    align.result        <- align( index          = base.string,
                                  readfile1      = read.1.files,
                                  readfile2      = read.2.files,
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
    return(rsubread.results)
}


get_power_analysis_result <- function() {

}

"
get_qc_final_obj()
setwd(output.bam.path)
save.image(image.data)
gene.counts <- get_rsubread_result()
save.image(image.data)
"









