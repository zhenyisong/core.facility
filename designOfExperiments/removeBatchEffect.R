#---
# Harman is a PCA and constrained optimisation 
# based technique that maximises the removal of 
# batch effects from datasets, with the constraint 
# that the probability of overcorrection (i.e. removing 
# genuine biological signal along with batch noise) is kept 
# to a fraction which is set by the end-user.
#---

# Bioconductor package: HarmanData, Harman
# Osmond-McLeod et al. 2013, Oytam et al. 2016

# @reference
# Batch effects : ComBat or removebatcheffects (limma package) ?
# https://www.biostars.org/p/266507/
# 
# 1: Nygaard V, Rødland EA, Hovig E. Methods that remove batch effects while
# retaining group differences may lead to exaggerated confidence in downstream
# analyses. Biostatistics. 2016 Jan;17(1):29-39.
#---

pkgs      <- c('Harman','HarmanData','limma','tidyverse','Biobase')
load.libs <- lapply(pkgs, require, character.only = TRUE)

data(OLF)

eset <- ExpressionSet( assayData   = as.matrix(olf.data),
                       phenoData   = AnnotatedDataFrame(olf.info))
# Plot principal components labeled by treatment
plotMDS( eset, labels = pData(eset)[, 2], gene.selection = 'common')
exprs(eset) <- removeBatchEffect(eset, batch = pData(eset)[,2])

# Plot principal components labeled by treatment
plotMDS(eset, labels = pData(eset)[,1], gene.selection = 'common')

# Plot principal components labeled by batch
plotMDS(eset, labels = pData(eset)[,2], gene.selection = 'common')