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

#---
# RSeQC annotation files
#---

rRNA_hg38_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/hg38_rRNA.bed'
house_keeping_genes_hg38_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/hg38.HouseKeepingGenes.bed'
basic_genes_gencode_hg38_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/hg38_GENCODE_v24_basic.bed'

rRNA_mm10_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/mm10_rRNA.bed'
house_keeping_genes_mm10_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/mm10.HouseKeepingGenes.bed'
basic_genes_gencode_mm10_RSeQC='/wa/zhenyisong/reference/annotation/RSeQC/mm10_GENCODE_VM11_basic.bed'


#---
# picard annotation files
#---

#---
# see more help on picard
# http://broadinstitute.github.io/picard/command-line-overview.html
#---

#---
# how to generate the refflat file
#
# for mm10
# please download the corresponding file from here
# http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/
# http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
# then I renamed to refFlat_mm10.txt
# gunzip refFlat.txt.gz
# mv refFlat.txt refFlat_mm10.txt
# mv refFlat.txt refFlat_hg38.txt
#---

#---
# how to generate mm10_ribosome_interval_list.txt"
#
# see: https://www.biostars.org/p/120145/
# see:http://seqanswers.com/forums/showthread.php?p=136425
# You can find the intervals using the UCSC Table browser. 
# For this, you go to 
# http://genome.ucsc.edu/cgi-bin/hgTables
# select mm10 version, mamalian
# and then set group:all tables, table:rmsk, 
# and filter to "repClass (does match) rRNA" 
# then output it as a GTF file.
# :: please set file name here, otherwise,
# :: the file will be displayed in browser
#--- 

REFFlAT_mm10_UCSC_picard='/wa/zhenyisong/reference/annotation/picard/refFlat_mm10.txt'
RIBO_INTERVAL_LIST_mm10_picard='/wa/zhenyisong/reference/annotation/picard/mm10_ribosome_interval_list.txt'

#---
# hisat2 index file location
# and other required annotation files
#---

hisat2_index_path='/home/zhenyisong/data/reference/index'





raw_data_path='/home/zhenyisong/data/results/chenlab/xiaoning/data'
working_dir='/home/zhenyisong/data/results/chenlab/xiaoning/qc.step/hisat2'
cd ${working_dir}


#---
# module 1
# aim -- generate bam file from raw data,
#     -- genome read mapping
#---


## if unpiared the data, -U parameter will be used
##shopt -s nullglob

read_1_files=($(find ${raw_data_path} -name '*_R1.fq.gz'))
read_2_files=($(find ${raw_data_path} -name '*_R2.fq.gz'))
file_num=${#read_1_files[@]}

echo 'now module I completed'
: << 'EOF'
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
    samtools sort -@ ${threads} -m 4G -O bam \
                  -o ${base}.sorted.bam ${base}.bam 
    # samtools fails, I do not know why?
    #samtools index ${base}.sorted.bam
    picard BuildBamIndex INPUT=${base}.sorted.bam
    
done

EOF
echo 'MODULE I end'
#---
# module 2
# aim -- use the fastqc module to analyze the raw data
#     
#---

mkdir fastqc.results
for (( i=0; i<${file_num}; i++ ));
do
    R1=${read_1_files[$i]}
    R2=${read_2_files[$i]}
    fastqc -o  fastqc.results -f fastq \
           -t ${threads} ${R1} ${R2}   
done



#---
# module 3
# aim -- using the picard procedure to infer the QC
#     -- rRNA percentage etc.
#---
STRANDNESS='NONE'

for (( i=0; i<${file_num}; i++ ));
do
    R1=${read_1_files[$i]}
    base=`basename $R1`
    base=${base%_R1.fq.gz}

    RANDOM_FILE_NAME=`date '+%m-%d-%H-%M-%s'`-`uuidgen -t`
    samtools view -H ${base}.sorted.bam -o ${RANDOM_FILE_NAME}
    cut -s -f 1,4,5,7,9  ${RIBO_INTERVAL_LIST_mm10_picard} >> ${RANDOM_FILE_NAME}
    picard CollectRnaSeqMetrics REF_FLAT=${REFFlAT_mm10_UCSC_picard} \
                                RIBOSOMAL_INTERVALS=${RANDOM_FILE_NAME} \
                                STRAND_SPECIFICITY=${STRANDNESS} \
                                CHART_OUTPUT=null \
                                METRIC_ACCUMULATION_LEVEL=ALL_READS \
                                INPUT=${base}.sorted.bam \
                                OUTPUT=${base}.picard  \
                                ASSUME_SORTED=true
    rm  ${RANDOM_FILE_NAME}
done

#---
# module 3
# aim -- using the RSeQC procedure to infer the strandness
#---
for (( i=0; i<${file_num}; i++ ));
do
    R1=${read_1_files[$i]}
    base=`basename $R1`
    base=${base%_R1.fq.gz}
    infer_experiment.py -r ${basic_genes_gencode_mm10_RSeQC} \
                        -i ${base}.sorted.bam
    read_distribution.py  -i ${base}.sorted.bam \
                          -r ${basic_genes_gencode_mm10_RSeQC}
    geneBody_coverage.py -r ${house_keeping_genes_mm10_RSeQC} \
                         -i ${base}.sorted.bam \  
                         -o ${base}.coverage.RSeQC
done

source deactivate macs2