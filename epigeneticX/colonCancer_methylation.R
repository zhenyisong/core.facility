# @author Yisong Zhen
# @since  2018-07-27
# @update 2018-08-10
#---


# @parents
#  colonCancer_methylation.sh  
# see more at Github:
#    genomicsclass/colonCancerWGBS
# some codes are from the assignment in this class
#---

library(bsseq)
library(bsseqData)
library(bumphunter)
library(tidyverse)
library(stringr)

colonCancer_data_path <- file.path('/wa/zhenyisong/cardiodata/SRP028600')
colonCancer_cov_files <- list.files( path = colonCancer_data_path, 
                                     pattern = 'bismark.cov.gz') %>%
                         paste0(colonCancer_data_path,'/',.)
sample_names          <- list.files( path = colonCancer_data_path, 
                                     pattern = 'bismark.cov.gz') %>%
                         str_replace('_1_bismark_bt2_pe.bismark.cov.gz', '')   

colonCancer_cov_BSseq <- read.bismark( files       = colonCancer_cov_files,
                                       sampleNames = sample_names,
                                       rmZeroCov   = TRUE,
                                       strandCollapse = FALSE,
                                       fileType = 'cov',
                                       mc.cores = 4,
                                       verbose  = TRUE,
                                       BACKEND  = NULL)
dim(colonCancer_cov_BSseq)
"
setwd(colonCancer_data_path)
save.image('colonCancer.Rdata')
"
"
cov = getCoverage(colonCancerWGBS,type  = 'Cov')
m   = getCoverage(colonCancerWGBS, type = 'M')
"
getCoverage(colonCancer_cov_BSseq, type = 'Cov' ) %>% head()
getCoverage(colonCancer_cov_BSseq, type = 'M' ) %>% head()

"
dat <- dummyData()
# Enable parallelization
require(doParallel)
registerDoParallel(cores = 2)
# Find bumps
bumps <- bumphunter(dat$mat, design=dat$design, chr=dat$chr, pos=dat$pos,
cluster=dat$cluster, coef=2, cutoff= 0.28, nullMethod="bootstrap",
smooth=TRUE, B=250, verbose=TRUE,
smoothFunction=loessByCluster)
"