#!/bin/bash

#---
# <book>
# Computational exome and Genome Anlysis
# I download the data using the wget method
#
#---

# qsub /wa/zhenyisong/sourcecode/GWAS/gwas.gatk4.sh
#---

source ~/.bash_profile
source ~/.bashrc
source activate biotools

#$ -N Yisong.GWAS
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/humangenetics/data
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -o /home/zhenyisong/data/humangenetics/data/job.log
#$ -e /home/zhenyisong/data/humangenetics/data/error.log
###$ -l h_vmem=16G

set -eo