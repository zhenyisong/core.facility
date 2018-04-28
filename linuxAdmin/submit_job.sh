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
#$ -N mRNAPythonQC
#$ -V
#$ -w e
#$ -wd /wa/zhenyisong/results/chenlab/songli/mRNAhumanYs
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log

#script_path='/wa/zhenyisong/sourcecode/core.facility/expressionNGS/songli_mRNA_chenlab.R'
#R CMD BATCH ${script_path}

#script_path='/wa/zhenyisong/sourcecode/core.facility/qualiControl/quality_control_industry.py'
#working_path='/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs/pythonQC'
#data_path='/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs'
#python ${script_path} -n '.clean.fq.gz' -g 'hg38' -l 'PE' \
#                      -s 'NONE' -w ${working_path} \
#                      -d ${data_path}

script_path='/wa/zhenyisong/sourcecode/core.facility/expressionNGS/songli_mRNA_chenlab.py'
python ${script_path}

source deactivate biotools