#!/bin/bash

# qsub /wa/zhenyisong/sourcecode/core.facility/linuxAdmin/aligner_index.sh

# please see more reference at this page
#
# http://gridscheduler.sourceforge.net/htmlman/manuals.html
# http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/
# sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html
# 
# http://talby.rcs.manchester.ac.uk/~ri/
# _linux_and_hpc_lib/sge_intro.html
#---

#---
# These two lines tell SGE to run the script using the BASH 
# command-line interpreter/shell ¡ª see below for details 
#---

#$ -S /bin/bash

#---
# The name of the job. The name should follow the  "name"
# definition  in sge_types(1).  Invalid job names will be
# denied at submit time.
#---

#$ -N Yisong.Aligner

#---
# The -V option to qsub states that the job should 
# have the same environment variables as the 
# shell executing qsub (recommended)
#---

#$ -V

#---
# Specifies a validation level applied to the job  to  be
# submitted  (qsub,  qlogin,  and  qsh)  or the specified
# queued job (qalter). 
#---

#$ -w e

#---
# Execute  the  job  from  the  directory  specified   in
# working_dir.
#---

#$ -wd /home/zhenyisong/data/reference/index

#---
# Mail at beginning/end/on suspension
#---

#$ -m ea

#---
# Send mail to these users
#---

#$ -M zhenyisong@gmail.com

#$ -j yes
#$ -o job.log
#$ -e error.log

###$ -l h_vmem=8G

#set -e
#set -u
#set -o pipefail

#---
# I asked wenke, and confirmed the cluster management system
# I think we can search the net using google
# "sge qsub script example'
#---
source ~/.bash_profile
source ~/.bashrc
source activate macs2
unset PYTHONPATH

# genome reference location
# tar xvzf igenomes/Mus_musculus_UCSC_mm9.tar.gz
# tar xvzf igenomes/Homo_sapiens_UCSC_hg19.tar.gz 1> /dev/null
#---
#mm10_UCSC_genome='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
#hg38_UCSC_genome='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
#rn06_UCSC_genome='/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa'
mm09_UCSC_genome='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa'
hg19_UCSC_genome='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'



# 
#---
# hisat2 index
# hisat2-build [options]* <reference_in> <ht2_base>
# https://ccb.jhu.edu/software/hisat2/manual.shtml
#---

#hisat2-build ${mm10_UCSC_genome} mm10
#hisat2-build ${hg38_UCSC_genome} hg38
#hisat2-build ${rn06_UCSC_genome} rn6
hisat2-build ${hg19_UCSC_genome} hg19
hisat2-build ${mm09_UCSC_genome} mm9



# rsubread index
#subread-buildindex -o mm10 ${mm10_UCSC_genome}
#subread-buildindex -o hg38 ${hg38_UCSC_genome}
#subread-buildindex -o rn6  ${rn06_UCSC_genome}
subread-buildindex -o hg19  ${hg19_UCSC_genome}
subread-buildindex -o mm9   ${mm09_UCSC_genome}



# bwa index
# 