# @author Yisong Zhen
# @since  2018-07-27
# @update 2018-07-27
#---


# @parents
#  colonCancer_methylation.sh  
# see more at Github:
#    genomicsclass/colonCancerWGBS
#---

library(bsseq)
library(bsseqData)
library(bumphunter)
library(tidyverse)
library(stringr)

colonCancer_data_path <- file.path('/wa/zhenyisong/cardiodata/SRP028600')
colonCancer_cov_files <- list.files( path = colonCaner_data_path, 
                                     pattern = 'bismark.cov.gz') %>%
                         paste0(colonCancer_data_path,'/',.)
sample_names          <- list.files( path = colonCaner_data_path, 
                                     pattern = 'bismark.cov.gz') %>%
                         str_replace('_1_bismark_bt2_pe.bismark.cov.gz', '')   
#setwd(colonCancer_data_path)
colonCancer_cov_BSseq <- read.bismark( files = colonCancer_cov_files,
                                       sampleNames = sample_names,
                                       rmZeroCov = FALSE,
                                       strandCollapse = FALSE,
                                       fileType = 'cov',
                                       mc.cores = 4,
                                       verbose  = TRUE,
                                       BACKEND  = NULL)
dim(colonCancer_cov_BSseq)
getCoverage(colonCancer_cov_BSseq) %>% head()