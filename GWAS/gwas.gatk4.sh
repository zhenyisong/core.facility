#!/bin/bash

#---
# <book>
# Computational exome and Genome Anlysis
# I download the data using the wget method
#
#---

# qsub /wa/zhenyisong/sourcecode/core.facility/GWAS/gwas.gatk4.sh
#---

source ~/.bash_profile
source ~/.bashrc
source activate macs2

#$ -S /bin/bash
#$ -N Yisong.GWAS
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/humangenetics/data
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
###$ -l h_vmem=16G

#--
# failed,
# I do not figure out why, it seemes to be related to shell enviroment
# at the HPC.
#---

#PS1='\[\e[1;32m\]HPC \[\e[1;33m\]\W [\@]\[\e[1;34m\] \[\e[1;35m\]$ \[\e[0m\]'
#
#set -e
#set -u 
#set -o pipefail

threads=3

#---
# data source I.
# Corpasome
# 1: Glusman G, Cariaso M, Jimenez R, Swan D, Greshake B, Bhak J, Logan DW, Corpas 
# M. Low budget analysis of Direct-To-Consumer genomic testing familial data.
# F1000Res. 2012 Jul 16;1:3.
# https://f1000research.com/articles/1-3/v1
# see methods section
# the book: p48.
#---

#---
# data source II.
# GIAB raw data
# Genome in a Bottle Consortium (GIAB) 
# https://github.com/genome-in-a-bottle
# the book: p49.
#---


#---
# igneomes
# hg38
#---
hg38='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
bwa_index='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome'
GIAB_raw_data_path='/wa/zhenyisong/humangenetics/data/Garvan_NA12878_HG001_HiSeq_Exome'
cd ${GIAB_raw_data_path}

# trim raw data
# /home/zhenyisong/data/miniconda3/envs/macs2/share/trimmomatic/adapters
#---

echo 'now initiat GWAS project'
: << 'EOF'
trimmomatic PE -threads ${threads} -trimlog NIST7035.log \
               NIST7035_TAAGGCGA_L001_R1_001.fastq.gz \
               NIST7035_TAAGGCGA_L001_R2_001.fastq.gz \
               NIST7035_trimmed_R1_paired.fastq.gz \
               NIST7035_trimmed_R1_unpaired.fastq.gz \
               NIST7035_trimmed_R2_paired.fastq.gz \
               NIST7035_trimmed_R2_unpaired.fastq.gz \
               ILLUMINACLIP:./adapters/NexteraPE-PE.fa:2:30:10 \
               LEADING:3 \
               TRAILING:3 \
               SLIDINGWINDOW:4:15 \
               MINLEN:36

EOF

echo 'finish first step'

#ln -s ${bwa_index}/* ./
read_group_info='@RG\tID:rg1\tSM:NA12878\tPL:illumina\tLB:lib1\tPU:H&AP8ADXX:1:TAAGGCGA'
bwa mem -t ${threads} -R '@RG\tID:rg1\tSM:NA12878\tPL:illumina\tLB:lib1\tPU:H&AP8ADXX:1:TAAGGCGA' \
           genome.fa \
           NIST7035_trimmed_R1_paired.fastq.gz \
           NIST7035_trimmed_R2_paired.fastq.gz  > NIST7035_aln.sam
samtools view -Sb NIST7035_aln.sam > NIST7035_aln.bam
picard ValidateSamFile INPUT=NIST7035_aln.bam MODE=SUMMARY
picard SortSam INPUT=NIST7035_aln.bam OUTPUT=NIST7035_sorted.bam SORT_ORDER=coordinate
           