#!/bin/bash
# @author Yisong Zhen
# @since  2018-06-26
# @update 2018-06-26
#---

# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/colonCancer_methylation.sh

#----
# HPC parameters for Sun Grid
#$ -S /bin/bash
#$ -N colonWGBSseq
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/cardiodata/SRP028600
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
###$ -l h_vmem=16G
#---


#---
# Please see the script at Github
# and HavardX course at Edx:
# HarvardX PH525.7x Data Analysis for Life Sciences 7 
# Case Studies in Functional Genomics
# Week-2_DNA_methylation
# colonCancerWGBS/scripts/createObject.Rmd 
#---
source ~/.bash_profile
source ~/.bashrc
source activate biotools
unset PYTHONPATH

hg38_igenome='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'

WGBS_index='/wa/zhenyisong/reference/WGBS/human'
raw_data='/wa/zhenyisong/cardiodata/SRP028600'
threads=2

if [ ! -d "${WGBS_index}" ]; then
    mkdir -p ${WGBS_index}
fi

cd ${WGBS_index}

if [ ! -f "genome.fa" ]
then
    ln -s ${hg38_igenome} ./
fi

bismark_genome_preparation --bowtie2 ./

cd ${raw_data}

find ./ -type f -name '*.sra' | xargs -n 1 -P 4 -I{} fastq-dump --gzip {}