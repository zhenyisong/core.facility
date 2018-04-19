#!/bin/bash

#---
#@author Yisong Zhen
#@since  2018-04-08
#update  2018-04-08
#---

# qsub /home/zhenyisong/data/sourcecode/core.facility/expressionNGS/songli_mRNA_QC.sh

#----
# HPC parameters for Sun Grid
#$ -N songlmRNAQC
#$ -S /bin/bash
#$ -w e
#$ -wd /home/zhenyisong/data/results/chenlab/songli/mRNAhumanYs
#$ -m ea
#$ -M zhenyisong@hotmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
#---


source ~/.bashrc
source ~/.bash_profile
source activate biotools
unset PYTHONPATH

PYTHON_QC='/home/zhenyisong/data/sourcecode/core.facility/qualiControl/quality_control_industry.py'

python ${PYTHON_QC} -n '.fq.gz' -g 'hg38' -l 'PE' -s 'NONE' -t 4

source deactivate biotools