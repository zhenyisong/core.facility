#!/bin/bash
# @author Yisong Zhen
# @since 2018-01-01
# this module is developed according to the requirement
# by the platform director Jingzhou Chen (Dr.) in the 
# mail corrrespondance.
#
# quality control (QC) pipeline to compelte or speed up the
# pre-QC process in B201 platform;
# most of the QC softwares were installed in macs2 enviroment
# I plan to use the linux shell to finish the job on QC-check
# already have several pipelines:
# snakemake?
# 
#

# I have to call the HPC parameters to use the HPC resource

# qsub /wa/zhenyisong/sourcecode/core.facility/qualiControl/quality_control.sh


#$ -N Yisong.QC
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/results/chenlab/xiaoning/qc.step/hisat2
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
###$ -l h_vmem=8G

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

# serveral choice to determine which kind of QC strategy to use
# use the switch phrase to pre-determine the setting.

#---
# define the computation power
# threads
#---

threads=3

# QC type
# mRNA-seq
# ChIP-seq
# lncRNA-seq
# miRNA-seq
# DNA-seq
#

# call upon the email messaging script to deliver the 
# final output is completed and the user should to notice
# the task is already finished. Please move on to the next steps.
#


#---
# annotation and genome files location
# location in linux platform is fixed
# the reference genome sets for various spieces
# were downloaded from iGenomes
#---
mm10_UCSC_genome='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
hg38_UCSC_genome='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
mm10_UCSC_GTF='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
hg38_UCSC_GTF='/home/zhenyisong/data/reference/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'

# mRNA-seq QC check list
# fasqc module
# RSeQC
# I download the required annotation files from 
# the author's website at
# http://rseqc.sourceforge.net/
# hg38 & mm10
# including the refseq (gencode, basic), house-keeping genes, and rRNA
# The annotation files were saved at specific folder created for RSeQC
#---

rRNA_hg38_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/hg38_rRNA.bed'
house_keeping_genes_hg38_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/hg38.HouseKeepingGenes.bed'
basic_genes_gencode_hg38_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/hg38_GENCODE_v24_basic.bed'

rRNA_mm10_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/mm10_rRNA.bed'
house_keeping_genes_mm10_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/mm10.HouseKeepingGenes.bed'
basic_genes_gencode_mm10_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/mm10_GENCODE_VM11_basic.bed'


#---
# hisat2 index file location
# and other required annotation files
#---

hisat2_index_path='/home/zhenyisong/data/reference/index'



#---
# module 1
# aim -- generate bam file from raw data,
#     -- genome read mapping
#---


## if unpiared the data, -U parameter will be used
##shopt -s nullglob

raw_data_path='/home/zhenyisong/data/results/chenlab/xiaoning/data'
working_dir='/home/zhenyisong/data/results/chenlab/xiaoning/qc.step/hisat2'
cd ${working_dir}
read_1_files=($(find ${raw_data_path} -name '*_R1.fq.gz'))
read_2_files=($(find ${raw_data_path} -name '*_R2.fq.gz'))
file_num=${#read_1_files[@]}


for (( i=0; i<${file_num}; i++ ));
do
    R1=${read_1_files[$i]}
    R2=${read_2_files[$i]}
    base=`basename $R1`
    base=${base%_R1.fq.gz}
    hisat2 -p ${threads} --dta --fr  \
           -x ${hisat2_index_path}/mm10 -1 ${R1} -2 ${R2} -S $base.sam
    samtools view -Sb -h -@ ${threads} -O BAM \
                  -T ${mm10_UCSC_genome} -o ${base}.bam  ${base}.sam  
done

source deactivate macs2
#---
# module 2
# aim -- using the picard procedure to infer the QC
#     -- rRNA percentage etc.
#---

#---
# module 3
# aim -- using the RSeQC procedure to infer the strandness
#---
