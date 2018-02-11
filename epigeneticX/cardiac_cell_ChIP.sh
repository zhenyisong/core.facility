#!/bin/bash


# qsub /home/zhenyisong/data/sourcecode/core.facility/epigeneticX/cardiac_cell_ChIP.sh


#----
# HPC parameters for Sun Grid
#$ -S /bin/bash
#$ -N Yisong.ChIPseq
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




raw_data_path='/wa/zhenyisong/cardiodata/GSE52386' 


#---
# script parameter setting
# output result dir and threads
#---

mapping_bowtie2_results='/wa/zhenyisong//cardiodata/GSE52386/bowtie2'
threads=6


#---
# index file paths
# and annotation files
#---
mm10_bowtie2_index='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index'


all_raw_data=($raw_data_path/*.fastq)
file_number=${#all_raw_data[@]}




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

#

if [ ! -d "${mapping_bowtie2_results}" ]; then
    mkdir -p "${mapping_bowtie2_results}"
    cd "${mapping_bowtie2_results}"
else
    cd "${mapping_bowtie2_results}"
fi


if [ ! -f 'genome.4.bt2' ]; then
    ln -s "${mm10_bowtie2_index}"/genome.* ./
fi
#---
# no need to reconstruct the BWA index
# ln -s $mm10 ./
# bwa index genome.fa > bwa.log 2>> bwa.log
#---

all_raw_data=($raw_data_path/*.fastq)
file_number=${#all_raw_data[@]}

echo 'first module'
: << 'EOF'
for (( i=0; i<$((file_number)); i++ ));
do
    filename=${all_raw_data[$i]}
    base=`basename ${filename}`
    base=${base%.fastq}
    bowtie2 -q -p ${threads} -x ${mm10_bowtie2_index}/genome -U $filename | \
    samtools view -bSh -@ ${threads} - | \
    picard SortSam INPUT='/dev/stdin' SORT_ORDER=coordinate OUTPUT="${base}.bam"
done


# merge two replicates in ChIP-seq data
#---

samtools merge -@ ${threads} -h SRR1029874.bam -O BAM SRR1029001.bam SRR1029874.bam SRR1029876.bam
samtools merge -@ ${threads} -h SRR1029875.bam -O BAM SRR1029002.bam SRR1029875.bam SRR1029877.bam

EOF

echo 'now complete the first step'


heart_whole_bams=($(find . -maxdepth 1 -type f -name "*.bam"))
index_array=($(seq 22 1 33))
index_array+=(48 49)
index_num=${#index_array[@]}

echo 'now the second step'
: << 'EOF'

for (( i=0; i<$((index_num)); i++ ));
do
    index=${index_array[$i]}
    bam_file=${heart_whole_bams[$index]}
    base=${bam_file%.bam}
    bamToBed -i ${bam_file} > ${base}.bed
done

EOF

echo 'now complete the second step'

for (( i=0; i<$((index_num)); i++ ));
do
    index=${index_array[$i]}
    treat=${heart_whole_bams[$index]}
    treat_base=${treat%.bam}
    i=$((i+1))
    index=${index_array[$i]}
    control=${heart_whole_bams[$index]}
    control_base=${control%.bam}
    epic --treatment "${treat_base}.bed"  --control "${control_base}.bed" --number-cores ${threads} \
         --genome mm10  --bed ${treat_base} > ${treat_base}.csv
done




