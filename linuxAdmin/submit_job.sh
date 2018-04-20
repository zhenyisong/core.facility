#!/bin/bash

# qsub /wa/zhenyisong/sourcecode/core.facility/linuxAdmin/submit_job.sh

#---
# I asked wenke, and confirmed the cluster management system
# I think we can search the net using google
# "sge qsub script example'
#---
source ~/.bash_profile
source ~/.bashrc
source activate biotools
unset PYTHONPATH


#$ -S /bin/bash
#$ -N Yisong.Xiaon
#$ -V
#$ -w e
#$ -wd /wa/zhenyisong/results/chenlab/songli/mRNAhumanYs
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log

script_path='/wa/zhenyisong/sourcecode/core.facility/expressionNGS/songli_mRNA_chenlab.R'
R CMD BATCH ${script_path}

source deactivate biotools