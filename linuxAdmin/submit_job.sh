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
#$ -wd /home/zhenyisong/data/sourcecode/ciona.network/code
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log

#script_path='/wa/zhenyisong/sourcecode/core.facility/expressionNGS/songli_phenotype_chenlab.R'
#R CMD BATCH ${script_path}

#script_path='/wa/zhenyisong/sourcecode/core.facility/qualiControl/code/quality_control_industry.py'
#working_path='/wa/zhenyisong/results/chenlab/songli/phenotype/pythonQC'
#data_path='/wa/zhenyisong/results/chenlab/songli/phenotype'
#python ${script_path} -n '.clean.fq.gz' -g 'hg38' -l 'PE' \
#                      -s 'NONE' -w ${working_path} \
#                      -d ${data_path}
#
#script_path='/wa/zhenyisong/sourcecode/core.facility/expressionNGS/songli_mRNA_chenlab.py'
#python ${script_path}

script_path='/wa/zhenyisong/sourcecode/ciona.network/code/ciona_working_code.py'
python ${script_path}

source deactivate biotools