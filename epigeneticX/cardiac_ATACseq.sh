#!/bin/bash
# @author  Yisong Zhen
# @since   2018-06-07
# @update  2018-06-07


# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/cardiac_ATACseq.sh


#----
# HPC parameters for Sun Grid
#$ -S /bin/bash
#$ -N Yisong.ATACseq
#$ -V
#$ -w e
#$ -wd /wa/zhenyisong/cardiodata/SRP101479
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
# GSE95763
# SRP101479
# I download the raw data using Aspera
# PMID: 28733351
# Quaife-Ryan GA, Sim CB, Ziemann M, Kaspi A et al. 
# Multicellular Transcriptional Analysis of Mammalian 
# Heart Regeneration. Circulation 2017 Sep 19;136(12):1123-1139.
# tissue source ; whole heart?
#--


threads=6

#---
# move data from subdir and decompress them.
# find ./ -type f -name '*.sra' | xargs -n 1 -P 4 -I{} fastq-dump --gzip {}
#---
raw_data_path='/wa/zhenyisong/cardiodata/SRP101479' 


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

mapping_bwa_results='/wa/zhenyisong/cardiodata/SRP101479/bwa'

mm10_bwa_index='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/BWAIndex'

all_raw_data=($raw_data_path/*.fastq.gz)


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
# bwa index -a bwtsw genome.fa > bwa.log 2>> bwa.log
# origianl protocol from the Cell paper
# S1, pdf= page12
# bwa aln -t 6 -l 25 mm9 sample:fastq:gz
#---



for filename in ${all_raw_data[@]}; 
do
    base=`basename ${filename}`
    base=${base%.fastq.gz}
    bwa mem -M -t 4 genome.fa ${raw_data_path}/${base}.fastq.gz > ${base}.bam
    picard SortSam INPUT=${base}.bam  OUTPUT=${base}.sorted.bam SORT_ORDER=coordinate
    picard MarkDuplicates INPUT=${base}.sorted.bam ASSUME_SORTED=true \
                          OUTPUT=${base}.dedup.sorted.bam METRICS_FILE=${base}.dedup.sorted.txt \
                          REMOVE_DUPLICATES=true 
done


for filename in ${all_raw_data[@]}; 
do
    base=`basename "${filename}"`
    base=${base%.fastq.gz}
    macs2 callpeak --treatment ${base}.dedup.sorted.bam \
          --format BAM --gsize mm --name ${base} \
          --keep-dup all --bdg --qvalue 0.01
done

source deactivate macs2
