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
#

# serveral choice to determine which kind of QC strategy to use
# use the switch phrase to pre-determine the setting.

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

source activate macs2

#---
# annotation and genome files location
# location in linux platform is fixed
# the reference genome sets for various spieces
# were downloaded from iGenomes
#---
mm10_UCSC_genome='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
mm9_UCSC_genome=''
hg19_UCSC_genome=''
hg38_UCSC_genome='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
mm10_UCSC_GTF='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
mm9_UCSC_GTF=''
hg19_UCSC_GTF=''
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


