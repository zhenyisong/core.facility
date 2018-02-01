﻿#!/bin/bash


# qsub /home/zhenyisong/data/sourcecode/core.facility/epigeneticX/cardiac_cell_ChIP.sh


#----
# HPC parameters for Sun Grid
#$ -S /bin/bash
#$ -N Yisong.MACS2
#$ -V
#$ -w e
#$ -wd /wa/zhenyisong/cardiodata/GSE52386
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
source activate macs2
unset PYTHONPATH


#---
# 
# the original data was downloaded from GEO
# GSE52386
# SRP/SRP033/SRP033009
# I download the raw data using Aspera
# PMID: 24360275 
# tissue source ; whole heart?
#--

#---
# using conda to install 
# conda install sra-tools
# vdb-config -i
# see this to edit the saving path
# https://www.biostars.org/p/159950/
# now download the raw data
# nohup prefetch SRR1611184&
#---

threads=6


raw_data_path='/wa/zhenyisong/cardiodata/GSE52386' 


#---
# munaully confirm if the bam file is ordered or not
# samtools view -H SRR1029856.bam
# the result should contain the specific FLAG
# SO:coordinate
# please see  more :
# https://www.biostars.org/p/5256/
#---
#---
# this path is linked to the BWA indexed file from
# igenome mm10 UCSC version
#---

mapping_bwa_results='/wa/zhenyisong/cardiodata/GSE52386/bwa'

mm10_bwa_index='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/BWAIndex'

all_raw_data=($raw_data_path/*.fastq)
file_number=${#all_raw_data[@]}


#

if [ ! -d "${mapping_bwa_results}" ]; then
    mkdir -p "${mapping_bwa_results}"
    cd "${mapping_bwa_results}"
else
    cd "${mapping_bwa_results}"
fi


if [ ! -f 'genome.fa.ann' ]; then
    ln -s "${mm10_bwa_index}"/genome.* ./
fi



#---
# how to construct the bwa index files
# bwa index genome.fa > bwa.log 2>> bwa.log
#---

all_raw_data=($raw_data_path/*.fastq)
file_number=${#all_raw_data[@]}


echo 'my step grid'
: << 'EOF'

for (( i=0; i<$((file_number)); i++ ));
do
    filename=${all_raw_data[$i]}
    base=`basename ${filename}`
    base=${base%.fastq}
    bwa aln -t ${threads} -l 25  ${mm10_bwa_index}/genome.fa $raw_data_path/${base}.fastq | \
    bwa samse ${mm10_bwa_index}/genome.fa - $raw_data_path/${base}.fastq | \
    samtools view -Shb -@ ${threads} - | \
    samtools sort -@ ${threads} -m 1G -o ${base}.bam -
done


#---
# merge two replicates in ChIP-seq data
# SRR1029001.bam & SRR1029002.bam
# are pseudo file names for next convenience
# and avoid to erase the alignment results.
# SRR1029001.bam will be the control
# SRR1029002.bam
#---

samtools merge -@ ${threads} -h SRR1029874.bam -O BAM SRR1029001.bam SRR1029874.bam SRR1029876.bam
samtools merge -@ ${threads} -h SRR1029875.bam -O BAM SRR1029002.bam SRR1029875.bam SRR1029877.bam

EOF

echo 'fish step chunck'


heart_whole_bams=($(find . -maxdepth 1 -type f -print0 -name "*.bam"))

bam_file_num=${#heart_whole_bams[@]}

for (( i=0; i<$((bam_file_num)); i++ ));
do
    treat=${heart_whole_bams[$i]}
    i=$((i+1))
    control=${heart_whole_bams[$i]}
    base=${treat%.bam}
    macs2 callpeak --treatment ${treat} --control ${control} \
          --format BAM --gsize mm --name ${base} --bdg --qvalue 0.01
done
