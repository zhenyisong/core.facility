#!/bin/bash

#---
# <book>
# Computational exome and Genome Anlysis
# I download the data using the wget method
#
#---

# qsub /wa/zhenyisong/sourcecode/core.facility/GWAS/gwas.gatk4.sh
#---

#----
# HPC parameters for Sun Grid
#$ -S /bin/bash
#$ -N Yisong.GATK4
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/humangenetics/data
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
###$ -l h_vmem=16G
#---

source ~/.bash_profile
source ~/.bashrc
source activate macs2
unset PYTHONPATH

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

threads=6

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
#
# data source II.
# GIAB raw data
# Genome in a Bottle Consortium (GIAB) 
# https://github.com/genome-in-a-bottle
# the book: p49.
# nohup wget -c -t 0 
# ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/
# Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz&
# nohup wget -c -t 0 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/
# Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz &
#
#---

#---
# data source III
# https://github.com/bahlolab/bioinfotools/blob/master/GATK/resource_bundle.md
# Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#---

Mills_standard_vcf='/wa/zhenyisong/reference/annotation/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
Mills_standard_vcf_index='/wa/zhenyisong/reference/annotation/GATK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
Mill_standard_vcf_filename='Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
Mills_standard_vcf_index_filename='Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'

#---
# igneomes
# hg38
#---
hg38='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
bwa_index_path='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex'
GIAB_raw_data_path='/wa/zhenyisong/humangenetics/data/Garvan_NA12878_HG001_HiSeq_Exome'
cd ${GIAB_raw_data_path}

# trim raw data
# /home/zhenyisong/data/miniconda3/envs/macs2/share/trimmomatic/adapters
#---

echo 'now first step'
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

if [ ! -f 'genome.fa.ann' ];  then
    ln -sf "${bwa_index_path}"/genome.* ./
fi



source deactivate macs2
source activate biotools
unset PYTHONPATH

##read_group_info="@RG\tID:rg1\tSM:NA12878\tPL:illumina\tLB:lib1\tPU:H&AP8ADXX:1:TAAGGCGA"
##bwa mem -t ${threads} -R ${read_group_info} \
##           genome.fa \
##           NIST7035_trimmed_R1_paired.fastq.gz \
##           NIST7035_trimmed_R2_paired.fastq.gz  > NIST7035_aln.sam
##samtools view -Shb -o NIST7035_aln.bam NIST7035_aln.sam
##gatk-launch ValidateSamFile --MODE SUMMARY --INPUT NIST7035_aln.bam 
##gatk-launch SortSam --INPUT NIST7035_aln.bam --OUTPUT NIST7035_sorted.bam --SORT_ORDER coordinate
##gatk-launch MarkDuplicates --INPUT NIST7035_sorted.bam --OUTPUT NIST7035_dedup.bam \
##                           --METRICS_FILE NIST7035.metrics --TAGGING_POLICY All
##gatk-launch BuildBamIndex --INPUT NIST7035_dedup.bam
##gatk-launch CreateSequenceDictionary --REFERENCE genome.fa --OUTPUT genome.dict
##samtools faidx genome.fa
##
##if [ ! -f "${Mill_standard_vcf_filename}"  -a  \
##     ! -f "${Mill_standard_vcf_index_filename}" ]; then
##      echo "file ${Mill_standard_vcf_filename} does not exists!"
##      ln -s ${Mills_standard_vcf} ./
##      ln -s ${Mills_standard_vcf_index} ./
##else
##    echo "now the Miller file is created using soft-link method!"
##fi
##
##
##
##
##gatk-launch HaplotypeCaller -R genome.fa -I NIST7035_dedup.bam \
##                            --genotyping-mode DISCOVERY -stand-call-conf 30 \
##                            -O raw_variants.vcf
##

##---
## hard filtering
##---

##filt_expr="QD < 2.0 || FS > 60.0 || MQ < 40.0 || \
##           MQRankSum < -12.5 || ReadPosRankSum < -8.0"
##gatk-launch  VariantFiltration -R genome.fa -V raw_variants.vcf \
##                               --filter-expression "${filt_expr}" \
##                               --filter-name 'snpFilter' \
##                               -O NA12878_filtered.vcf

# result raw eyes checking

#bcftools query -f '%FILTER\n' NA12878_filtered.vcf | wc -l
#
#bcftools query -f '%FILTER\n' NA12878_filtered.vcf | grep PASS | wc -l
# 

# the hg38 ncbi ,cannot download the files
# I give up the try and instead using ucsc to 
# jannovar download -d hg38/refseq
# this commnand canot be carried out in cluster node
# as node has no right to ftp the outside.
#  jannovar/manual/datasource.rst
#  jannovar/jannovar-cli/src/main/resources/default_sources.ini
#---

#jannovar download --data-source-list ~/data/sourcecode/core.facility/GWAS/hg38.ini -d hg38/refseq
# we have to start sslocal to download and generate the files from NCBI or EBI
# I manualyy downloaded the required files from the corresponding
# website used in 
# https://gist.github.com/holtgrewe/5517986d90fe551a76b6091340eda79e#file-patched_sources-ini-L176
# and saved to the corresponding dir.
#----

#---
# for the memory increase strategy
# https://stackoverflow.com/questions/29577726/increasing-java-heap-size-pace
#---
JAVA_TOOL_OPTIONS="-Xms5g -Xmx5g -XX:+UseConcMarkSweepGC -XX:-UseGCOverheadLimit -XX:MetaspaceSize=4g"
#---
# using the NCBI annotation version failed.
# 
#JAVA_TOOL_OPTIONS="-Xms2G -Xmx4G"
#jannovar download -d hg38/refseq  --http-proxy http://127.0.0.1:1080 \
#        --ftp-proxy http://127.0.0.1:1080
#jannovar download  -d hg38/refseq
#jannovar annotate-vcf -d data/hg38_refseq.ser -i NA12878_filtered.vcf \
#                     -o NA12878_annotated.vcf
# the above code is broken due to the memory limit in this small java program
# I give up any further try and used the UCSC version instead.
#---
jannovar download  -d hg38/ucsc
jannovar annotate-vcf -d data/hg38_ucsc.ser -i NA12878_filtered.vcf \
                      -o NA12878_annotated.vcf