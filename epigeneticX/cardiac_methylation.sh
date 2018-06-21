#!/bin/bash
# @author Yisong Zhen
# @since  2018-06-11
# @update 2018-06-19
#---

# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/cardiac_methylation.sh

#----
# HPC parameters for Sun Grid
#$ -S /bin/bash
#$ -N WGBSseq
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/cardiodata/SRP017503
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
###$ -l h_vmem=16G
#---

set -eo
#---
# I asked wenke, and confirmed the cluster management system
# I think we can search the net using google
# "sge qsub script example'
#---
source ~/.bash_profile
source ~/.bashrc
source activate biotools
unset PYTHONPATH


WGBS_index='/wa/zhenyisong/reference/WGBS/mouse'
raw_data='/wa/zhenyisong/cardiodata/SRP017503'
threads=2
#---
# conda install bismark
# cd ${WGBS_index}
# bismark_genome_preparation --bowtie2 ./
#---

#---
# The current Illumina protocol for BS-Seq is directional, 
# in which case the strands complementary
# to the original strands are merely theoretical 
# and should not exist in reality.
#---
cd ${raw_data}

#find ./ -type f -name '*.sra' | xargs -n 1 -P $threads -I{} fastq-dump --gzip {}

for filename in *.fastq.gz; 
do 
    bismark --multicore $threads --bowtie2 --bam ${WGBS_index} $filename
done