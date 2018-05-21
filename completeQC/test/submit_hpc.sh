#!/bin/bash

#---
# aims
#    this Snippet code is the bundled script for
#    the QC python pipeline to be submitted to
#    HPC at the Clinic Test Center.
# @author Yisong Zhen
# @since  2018-05-03
# @update 2018-05-16
#---


# qsub /wa/zhenyisong/sourcecode/core.facility/completeQC/test/submit_hpc.sh


#---
# I asked wenke, and confirmed the cluster management system
# I think we can search the net using google
# "sge qsub script example'
#---
source ~/.bash_profile
source ~/.bashrc
source activate biotools

#---
# the command line specifications
# for 
#$ -S /bin/bash
#$ -N MyQCJob
#$ -V
#$ -w e
#$ -wd /wa/zhenyisong/results/chenlab/songli/phenotype
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log

#---
# the script installatio path
#---

script_path='/home/zhenyisong/data/sourcecode/core.facility/expressionNGS/songli_phenotype_chenlab.R'

#---
# @parameter
#     working_path
# you need to specify the working path which the QC results will be saved 
# there.
#---
working_path='/home/zhenyisong/data/results/chenlab/songli/phenotype/pythonQC'

#---
# @parameter
#     data_path
# you need to specify the data path where the raw reads sequencing data
# will be called from.
#---
data_path='/home/zhenyisong/data/results/chenlab/songli/phenotype'
quality_control_industry.py -n '.clean.fq.gz' -g 'hg38' -l 'PE' \
                            -s 'NONE' -w ${working_path} \
                            -d ${data_path}

#R CMD BATCH ${script_path}

source deactivate biotools