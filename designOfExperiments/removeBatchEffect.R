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

pkgs <- c('Harman','HarmanData')
load.libs <- apply(pkgs, 1, )