#!/bin/bash

# qsub /wa/zhenyisong/sourcecode/core.facility/expressionNGS/xiaon_cluster.sh


#---
# I asked wenke, and confirmed the cluster management system
# I think we can search the net using google
# "sge qsub script example'
#---
source ~/.bash_profile
source ~/.bashrc
source activate biotools

#$ -S /bin/bash
#$ -N Yisong.Xiaon
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/results/chenlab
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
##$ -l h_vmem=8G

source_code_path='/home/zhenyisong/data/sourcecode/core.facility/expressionNGS'

#Rscript ${source_code_path}/xiaoning_chenlab.R
R CMD BATCH ${source_code_path}/xiaon_speed_chenlab.R

source deactivate biotools