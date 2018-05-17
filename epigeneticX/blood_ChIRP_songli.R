# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-05-17
# @parent  blood_ChIRP_songli.sh
#---

pkgs              <- c( 'tidyverse', 'GenomicRanges',
                        'GenomicAlignments', 'BiocParallel',
                        'Rsamtools','magrittr', 'DESeq2',
                        'stringr', 'PWMEnrich.Hsapiens.background',
                        'PWMEnrich',
                        'TxDb.Hsapiens.UCSC.hg38.knownGene',
                        'Homo.sapiens', 'ggbio', 'ChIPpeakAnno',
                        'org.Hs.eg.db', 'clusterProfiler',
                        'edgeR','DESeq2', 'openxlsx',
                        'BSgenome.Hsapiens.UCSC.hg38')
             
load.lib          <- lapply(pkgs, require, character.only = TRUE)

macs2.ChIRP.path          <- file.path('/wa/zhenyisong/results/chenlab/songli/ChIPRbwa')
macs14.bowtie2.ChIRP.path <- file.path('/home/zhenyisong/data/results/chenlab/songli/bowtie2')

working.env <- 'linux'
linux.path  <- macs2.ChIRP.path 
window.path <- file.path('D:\\yisong.data')
image.data  <- 'songliChIRP.Rdata'


# rm(list=ls())

switch( working.env, linux  = { setwd(linux.path);
                                load( file = image.data)},
                     window = { setwd(window.path);
                                load( file = image.data)} )
#---
# https://support.bioconductor.org/p/83599/
# to merge the GRange object?
#---

setwd(macs2.ChIRP.path)

z.ChIRP.hg38.df        <- read.delim( 'blood_z_peaks.xls', 
                                       comment.char = '#')
summary(z.ChIRP.hg38.df)

nova.ChIRP.hg38.df     <- read.delim( 'blood_nova_peaks.xls', 
                                       comment.char = '#')
summary(nova.ChIRP.hg38.df)

chromosome.set       <- paste('chr',c(1:22,'X','Y'), sep = '')  
nova.ChIRP.hg38.macs2  <- nova.ChIRP.hg38.df              %$%
                          { GRanges(  
                            seqname      = paste('chr',chr, sep = ''),
                            ranges       = IRanges( start = start, 
                                                    end   = end),
                            peak.length  = length,
                            abs.summit   = abs_summit,
                            pileup       = pileup,
                            log.pvalue   = X.log10.pvalue.,
                            foldChange   = fold_enrichment,
                            logqvalue    = X.log10.qvalue.,
                            macs2.name   = name ) } %>%
                          {.[seqnames(.) %in% chromosome.set]}
z.ChIRP.hg38.macs2     <- z.ChIRP.hg38.df                %$%
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
                            macs2.name   = name ) } %>%
                           {.[seqnames(.) %in% chromosome.set]}



two.studies <- suppressWarnings( findOverlaps( query   = nova.ChIRP.hg38.macs2, 
                                               subject = z.ChIRP.hg38.macs2, 
                                               type    = 'any') )
               
z.ChIRP.macs2.summit <- read.delim( 'blood_z_summits.bed', header = FALSE) %$%
                        { GRanges(  
                            seqname      = V1,
                            ranges       = IRanges( start = V2, 
                                                    end   = V3),
                            logqvalue    = V5,
                            macs2.name   = V4 ) } 
end(z.ChIRP.macs2.summit)   <- end(z.ChIRP.macs2.summit) + 1000
start(z.ChIRP.macs2.summit) <- start(z.ChIRP.macs2.summit) + 1000


bwa.z.bams                  <- c( 'ODD_clean.bam',
                                  'Even_clean.bam', 
                                  'blood.full.bam',
                                  'Input_clean.bam') %>%
                               Rsamtools::BamFileList(
                                         yieldSize = 7500000)


register( MulticoreParam( workers = 6) )
bwa.read.counts            <- summarizeOverlaps(  
                                  features      = z.ChIRP.hg38.macs2,
                                  reads         = bwa.z.bams, 
                                  ignore.strand = TRUE, 
                                  singleEnd     = FALSE ) %>% assay()

# multiBamSummary bins --bamfiles ODD_clean.bam Even_clean.bam -o results.temp.npz
# plotCorrelation --corData results.temp.npz --corMethod pearson \
#                   --outFileCorMatrix results.temp.txt --plotFile results.temp.pdf \
#                  --whatToPlot scatterplot 
       
ChIRP.hg38.annot     <- toGRanges( TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                   feature = 'gene')

bwa.normed.oe.counts <- bwa.read.counts[,1:2] %>%
                        cpm(normalized.lib.sizes = TRUE, log = TRUE)
    
bwa.counts.sd        <- apply(bwa.normed.oe.counts, 1, sd)
summary(bwa.counts.sd)

seqlevelsStyle(z.ChIRP.hg38.macs2) <- seqlevelsStyle(ChIRP.hg38.annot)

# do annotation by nearest TSS
blood.peaks.annot        <- z.ChIRP.hg38.macs2[ 
                                bwa.counts.sd < 0.6 &
                                mcols(z.ChIRP.hg38.macs2)$foldChange >=3 ] %>%
                            annotatePeakInBatch( ,
                                   AnnotationData = ChIRP.hg38.annot,
                                   select         = 'first')
       

peaks.hg38.gene.symbols  <- mapIds( org.Hs.eg.db, 
                                    keys      = blood.peaks.annot %>%
                                                mcols() %$% feature %>%
                                                unique(),
                                    column    = 'SYMBOL', 
                                    keytype   = 'ENTREZID', 
                                    multiVals = 'first') 

# this need to save the processed data in linux evn.
"
setwd(macs2.ChIRP.path)
save.image('songliChIRP.Rdata')
"

songli.ChIRP.table       <- peaks.hg38.gene.symbols %>% 
                            names() %>% unique() %>%
                            enrichKEGG( organism = 'human', 
                                        pvalueCutoff  = 0.05, 
                                        pAdjustMethod = 'none',
                                        qvalueCutoff  = 1) %>%
                            summary()
songli.ChIRP.pathways    <- songli.ChIRP.table %$%
                            {data.frame( kegg.pvalue  = -log(pvalue),
                                         kegg.pathway = Description )} %>%
                            ggplot( aes( x = reorder(kegg.pathway, kegg.pvalue), 
                                         y = kegg.pvalue))   + 
                            geom_bar( stat     = 'identity', 
                                      width    = 0.4, 
                                      position = position_dodge( width = 0.1), 
                                      size     = 20) +
                            theme( axis.text.x = element_text( angle = 60, 
                                                               hjust = 1, 
                                                               size  = 8 ),
                                   axis.text.y = element_text( hjust = 1, 
                                                               size  = 10) ) +
                            ylab('-log(pvalue)') + 
                            xlab('ChIRP.KEGG.ggplot') + 
                            coord_flip()
(songli.ChIRP.pathways)



# MF seleciton criteria
# peaks.hg38.gene.symbols %>% 
#                            names() %>% unique() %>% 
#                   enrichGO( OrgDb         = org.Hs.eg.db,
#                             ont           = 'MF',
#                             pAdjustMethod = 'none',
#                             pvalueCutoff  = 0.01,
#                             qvalueCutoff  = 1) %>%
#                   summary() %$%
#                   {data.frame( Golabels  = paste(ID, Description, sep = ' '),
#                                LogPvalue = - log(pvalue) ) } [c(1:9,12:16),]%>%

songli.ChIRP.MF.table <- peaks.hg38.gene.symbols %>% 
                            names() %>% unique() %>% 
                         enrichGO( OrgDb         = org.Hs.eg.db,
                                   ont           = 'MF',
                                   pAdjustMethod = 'none',
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 1) %>%
                         summary()
songli.ChIRP.BP.table <- peaks.hg38.gene.symbols %>% 
                            names() %>% unique() %>% 
                         enrichGO( OrgDb         = org.Hs.eg.db,
                                   ont           = 'BP',
                                   pAdjustMethod = 'none',
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 1) %>%
                         summary() 
songli.ChIRP.CC.table <- peaks.hg38.gene.symbols %>% 
                            names() %>% unique() %>% 
                         enrichGO( OrgDb         = org.Hs.eg.db,
                                   ont           = 'CC',
                                   pAdjustMethod = 'none',
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 1) %>%
                         summary()                   
songli.ChIRP.GO       <- songli.ChIRP.CC.table %$%
                         {data.frame( Golabels  = paste(ID, Description, sep = ' '),
                                      LogPvalue = - log(pvalue) ) } [c(1:9,12:16),]%>%
                         ggplot( aes( x = reorder(Golabels, LogPvalue), 
                                                    y = LogPvalue)) + 
                         geom_bar( stat = 'identity', width = 0.4, 
                                  position = position_dodge(width = 0.1), size = 20) +
                         theme( axis.text.x = element_text(angle = 60,hjust = 1, size = 8),
                                axis.text.y = element_text(hjust = 1, size = 10)) +
                         ylab('-log(pvalue)') + 
                         xlab('ChIRP.MF.ggplot2') + 
                         coord_flip()

(songli.ChIRP.GO)


#---
# the following is response to the songli's request
# to get the fine genome coordinate of these genes 
# from her list
# 2018-04-18
#---
macs2.z.annot.df         <- annotatePeakInBatch(z.ChIRP.hg38.macs2 ,
                                 AnnotationData = ChIRP.hg38.annot,
                                 select         = 'first') %>%
                            as.data.frame %>% cbind(bwa.read.counts)
    
                         
whole.peaks.gene.symbols <- mapIds( org.Hs.eg.db, 
                                    keys      = macs2.z.annot.df %$%
                                                feature,
                                    column    = 'SYMBOL', 
                                    keytype   = 'ENTREZID', 
                                    multiVals = 'first') 
whole.peaks.annotation   <- macs2.z.annot.df %>% cbind(whole.peaks.gene.symbols)

songli.selections        <- c( 'CFAP57', 'HORMAD1', 'NID1',
                               'NLRP3', 'SFTPD', 'ST14',
                               'LINC01220','TSKS','DUSP2',
                               'ACOXL', 'TMEM191B','ACSL1',
                               'EGR1','DDX43', 'LINC00689',
                               'SLC18A1','LZTS1','GLDC','CNTNAP3'
                              )
songli.df                <- whole.peaks.annotation %>% 
                            filter(whole.peaks.gene.symbols %in% songli.selections)
songli.df$whole %>% as.character() %>% unique()


#---
# find the enriched motifs
#---

data(PWMLogn.hg19.MotifDb.Hsap)

get_enriched_results <- function(cutoff) {
    enriched_report  <- whole.peaks.annotation %>%
                        filter(blood.full.bam >= cutoff) %$%
                        { GRanges(  
                          seqname      = seqnames,
                          ranges       = IRanges( start = start, 
                                                  end   = end),
                          strand       = '+',
                          peak.length  = peak.length,
                          abs.summit   = abs.summit,
                          pileup       = pileup,
                          log.pvalue   = log.pvalue,
                          foldChange   = foldChange,
                          logqvalue    = logqvalue,
                          macs2.name   = macs2.name ) } %>%
                        getSeq(x = BSgenome.Hsapiens.UCSC.hg38, names = .) %>%
                        motifEnrichment(PWMLogn.hg19.MotifDb.Hsap) %>%
                        groupReport()  
    return(enriched_report)

}

cutoff      <- seq(100, 400, 50)

report_list <- map(cutoff, get_enriched_results)

"
setwd(macs2.ChIRP.path)
save.image('songliChIRP.Rdata')
q('no')
"

                              

#setwd(window.path)
songli.ChIRP.wb <- createWorkbook()

#---
# if this in the linux env,
# the following code must be carried out
# otherwise, openslxs will be failed.
#---

"
zip.path <- '/wa/zhenyisong/miniconda3/envs/biotools/bin/zip'
Sys.setenv('R_ZIPCMD' = zip.path)
"




addWorksheet(songli.ChIRP.wb, 'KEGG pathway')
addWorksheet(songli.ChIRP.wb, 'GO-MF')
addWorksheet(songli.ChIRP.wb, 'GO-CC')
addWorksheet(songli.ChIRP.wb, 'GO-BP')
addWorksheet(songli.ChIRP.wb, 'Gene-symbol')
addWorksheet(songli.ChIRP.wb, 'songliChoice')
addWorksheet(songli.ChIRP.wb, 'wholePeaksAnnotation')


writeData(songli.ChIRP.wb, sheet = 1, songli.ChIRP.table )
writeData(songli.ChIRP.wb, sheet = 2, songli.ChIRP.MF.table )
writeData(songli.ChIRP.wb, sheet = 3, songli.ChIRP.CC.table )
writeData(songli.ChIRP.wb, sheet = 4, songli.ChIRP.BP.table )
writeData(songli.ChIRP.wb, sheet = 5, peaks.hg38.gene.symbols )
writeData(songli.ChIRP.wb, sheet = 6, songli.df )
writeData(songli.ChIRP.wb, sheet = 7, whole.peaks.annotation )
saveWorkbook(songli.ChIRP.wb, 'songliChIRP.2018-04-20.xlsx', overwrite = TRUE)