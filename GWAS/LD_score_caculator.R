# @author Yisong Zhen
# @since   2018-10-19
# @reference
#     Edx:
#     KyotoUx: 005x
#     Introduction to Statistical Methods for Gene Mapping
#     books:
#     Applied Statistical Genetics with R: For Population-based Association Studies
#     Page. 65
#---

library(tidyverse)
library(genetics)
human_data.path <- file.path('/wa/zhenyisong/sourcecode/GWAS/data')
setwd(human_data.path)
fms   <- read.delim('FMS_data.txt', header = T, sep = '\t')

# throw out errors.
"
raw.data <- read_delim('FMS_data.txt', delim = '\t',
                        col_names = TRUE, comment = '#')
"