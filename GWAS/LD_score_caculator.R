# @author Yisong Zhen
# @since   2018-10-22
# @reference
#     Edx:
#     KyotoUx: 005x
#     Introduction to Statistical Methods for Gene Mapping
#     books:
#     Applied Statistical Genetics with R: For Population-based Association Studies
#     Page. 65
#     http://www.stat-gen.org/
#---

library(tidyverse)
library(genetics)
human_data.path <- file.path('/wa/zhenyisong/sourcecode/core.facility/GWAS/data')
setwd(human_data.path)
fms   <- read.delim('FMS_data.txt', header = T, sep = '\t')
attach(fms)
Actn3Snp1 <- genotype(actn3_r577x,    sep = '')
Actn3Snp2 <- genotype(actn3_rs540874, sep = '')
LD(Actn3Snp1, Actn3Snp2)$"D'"

# throw out errors.
"
raw.data <- read_delim('FMS_data.txt', delim = '\t',
                        col_names = TRUE, comment = '#')
"

hgdp <- read.delim('HGDP_AKT1.txt', header = T, sep = '\t')
attach(hgdp)
Akt1Snp1 <- AKT1.C0756A
ObsCount <- table(Akt1Snp1)

ObsDat <- matrix(c(48,291,291,724), byrow = T, 2, 2)
chisq.test(ObsDat, correct = FALSE)