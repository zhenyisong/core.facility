#---
# @author Yisong Zhen
# @since  2017-12-29
# @update 2017-12-29
# cheat sheet
# https://conda.io/docs/user-guide/cheatsheet.html
# 
# please visit the website at
# https://bioconda.github.io/index.html
# how to install bioconda
# Bioconda is a channel
# or see the recent paper about Bioconda project
#
# 'Bioconda: A sustainable and comprehensive software distribution
#  for the life sciences'
#---


#---
# https://conda.io/docs/user-guide/tutorials/build-pkgs-skeleton.html
# Building conda packages with conda skeleton

# https://conda.io/docs/user-guide/index.html
# conda user guide

# I first downloaded the linux version 
# https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# this is full version of conda. > 500M
# this version of conda was discared in my current settting
#---

#---
# how to uninstall
# Uninstalling Anaconda or Miniconda
#
# remove previous legacy
# rm -rf ~/.condarc ~/.conda ~/.continuum
# and remove whole dir /wa/zhenyisong/anaconda3
#---

#https://conda.io/miniconda.html
bash Miniconda3-latest-Linux-x86_64.sh
# yes & /wa/zhenyisong/miniconda3

#---
# this is default, no need to be specified
# conda config --add channels defaults
# I added two following channels for
# bioinforamtics task
#---
conda config --add channels conda-forge
conda config --add channels bioconda

conda create --name biotools
source activate biotools
conda env list

#---
# no need
# Sys.setenv(CURL_CA_BUNDLE = '/wa/zhenyisong/anaconda3/envs/biotools/ssl/cacert.pem')
#
# conda search r-essentials
#---

conda install r-essentials
conda install curl=7.49.1
conda install bioconductor-ggbio bioconductor-gviz

conda create --name macs2
source activate macs2
conda install python=2.7 bwa samtools macs2

# source activate biotools
# now run R
source activate macs2
conda list --revisions > macs2.conda.history
conda list --explicit > macs2.conda.env
conda env export > macs2.conda.yml