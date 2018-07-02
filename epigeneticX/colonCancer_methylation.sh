#!/bin/bash
# @author Yisong Zhen
# @since  2018-06-26
# @update 2018-07-02
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
source ~/.bashrc
source ~/.bash_profile
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

#---
# Building the bisulfite genome indexes using Bowtie 2
#---
#bismark_genome_preparation --bowtie2 ./


#---
# decompress the SRA files in the folder
#---
cd ${raw_data}

#find ./ -type f -name '*.sra' | xargs -n 1 -P ${threads} -I{} fastq-dump --split-files --gzip {}

#---
# Using the Bismark read aligner
#---

all_raw_data_read_1=($raw_data/*_1.fastq.gz)
all_raw_data_read_2=($raw_data/*_2.fastq.gz)
file_number=${#all_raw_data_read_1[@]}




for (( i=0; i<$((file_number)); i++ ));
do
    filename=${all_raw_data_read_1[$i]}
    base=`basename ${filename}`
    base=${base%_1.fastq.gz}
    bismark --multicore ${threads} --bowtie2 --bam ${WGBS_index} \
            -1 ${base}_1.fastq.gz -2 ${base}_2.fastq.gz
done

#---
# The --bedGraph parameter generates coverage files which 
# are useful as input for the R/Bioconductor packge bsseq later on. 
# In addition, an M-bias report will be generated and an overall count report.
#---

#for filename in *.bam;
#do
#    bismark_methylation_extractor -p --no_overlap --comprehensive \ 
#                                  --multicore ${threads} --buffer_size 5G --bedGraph --counts --gzip \
#                                  ${filename}SRR949215_1_val_1.fq_bismark_bt2_pe.bam
#done
