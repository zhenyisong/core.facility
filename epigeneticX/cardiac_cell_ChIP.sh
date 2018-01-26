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


for (( i=0; i<$((file_number)); i++ ));
do
    filename=${all_raw_data[$i]}
    base=`basename ${filename}`
    base=${base%.fastq}
    bowtie2 -q -p ${threads} -x ${mm10_bowtie2_index}/genome -U $filename | \
    samtools view -Shb -@ ${threads} - | \
    picard SortSam INPUT='/dev/stdin' OUTPUT="${base}.bam"
done



exit 0

# merge two replicates in ChIP-seq data
#---

samtools merge -@ ${threads} -h SRR1029874.sam heart_11.5.treat.sorted.bam  SRR1029874.sorted.bam SRR1029876.sorted.bam
samtools merge -@ ${threads} -h SRR1029875.sam heart_11.5.control.sorted.bam SRR1029875.sorted.bam  SRR1029877.sorted.bam



mv SRR1029874.sorted.bam SRR1029874.sorted.bam.old
mv SRR1029876.sorted.bam SRR1029876.sorted.bam.old
mv SRR1029875.sorted.bam SRR1029875.sorted.bam.old
mv SRR1029877.sorted.bam SRR1029877.sorted.bam.old


# the epic
##new_bam_files=($mapping_results/*.sorted.bam)
##bam_file_num=${#new_bam_files[@]}
##for (( i=0; i<$((bam_file_num)); i++ ));
##do
##    filename=${new_bam_files[$i]}
##    base=`basename "${filename}"`
##    base=${base%.sorted.bam}
##    bamToBed -i ${base}.sorted.bam > ${base}.bed
##done
##
##heart_treat_beds=(heart_11.5.treat.bed  SRR1029878.bed SRR1029880.bed 
##                  SRR1029882.bed SRR1029884.bed SRR1029886.bed SRR1029888.bed)
##heart_control_beds=(heart_11.5.control.bed SRR1029879.bed SRR1029881.bed 
##                    SRR1029883.bed SRR1029885.bed SRR1029887.bed SRR1029889.bed)
##
##bed_file_num=${#heart_treat_beds[@]}
##
##for (( i=0; i<$((bed_file_num)); i++ ));
##do
##    treat_bed=${heart_treat_beds[$i]}
##    control_bed=${heart_control_beds[$i]}
##    base=${treat_bed%.bed}
##    epic --treatment ${treat_bed}  --control ${control_bed} --number-cores ${threads} \
##         --genome mm10  --bed ${base}_bed > ${base}.csv
##done
##



exit 0

heart_treat_bams=(heart_11.5.treat.sorted.bam  SRR1029878.sorted.bam SRR1029880.sorted.bam 
                  SRR1029882.sorted.bam SRR1029884.sorted.bam SRR1029886.sorted.bam SRR1029888.sorted.bam)
heart_control_bams=(heart_11.5.control.sorted.bam SRR1029879.sorted.bam SRR1029881.sorted.bam 
                    SRR1029883.sorted.bam SRR1029885.sorted.bam SRR1029887.sorted.bam SRR1029889.sorted.bam)

bam_file_num=${#heart_treat_bams[@]}

for (( i=0; i<$((bam_file_num)); i++ ));
do
    treat_bam=${heart_treat_bams[$i]}
    control_bam=${heart_control_bams[$i]}
    base=${treat_bed%.sorted.bam}
    macs2 callpeak --treatment ${treat_bam} --control ${control_bam} \
          --format BAM --gsize mm --name ${base} --bdg --qvalue 0.01
done
