#!/bin/bash
# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-03-22
#---


# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/blood_CHiRP_songli.sh


#----
# HPC parameters for Sun Grid
#$ -N Yisong.MACS2
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/results/chenlab/songli
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
###$ -l h_vmem=16G
#---

#---
# I asked wenke, and confirmed the cluster management system
# I think we can search the net using google
# "sge qsub script example'
# follwoing script is runing on the macs2 env by conda management
#
#---
source ~/.bash_profile
source activate macs2
unset PYTHONPATH


#---
# 
# the original data were uploaded from songli
# she has a project on CHiRP.
# I have already checked the data fingerprints
# using the command
# md5sum -c CleanFq_md5sum.TXT
# and all data are intact and in good state
#--


threads=6


raw_data_path='/wa/zhenyisong/results/chenlab/songli/CleanFq' 


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

mapping_bwa_results='/wa/zhenyisong/results/chenlab/songli/bwa'

hg38_bwa_index='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex'


#

if [ ! -d "${mapping_bwa_results}" ]; then
    mkdir -p "${mapping_bwa_results}"
    cd "${mapping_bwa_results}"
else
    cd "${mapping_bwa_results}"
fi


if [ ! -f 'genome.fa.ann' ]; then
    ln -s "${hg38_bwa_index}"/genome.* ./
fi



#---
# how to construct the bwa index files
# bwa index -a bwtsw genome.fa > bwa.log 2>> bwa.log
# origianl protocol from the Cell paper
# S1, pdf= page12
# bwa aln -t 6 -l 25 mm9 sample:fastq:gz
#---

all_raw_data=($raw_data_path/*.fq.gz)
file_number=${#all_raw_data[@]}


#echo 'my step grid'
#: << 'EOF'

for (( i=0; i<$((file_number)); i++ ));
do
    filename=${all_raw_data[$i]}
    base=`basename ${filename}`
    base=${base%.fq.gz}
    bwa aln -t ${threads} -l 25  ${hg38_bwa_index}/genome.fa $raw_data_path/${base}.fq.gz | \
    bwa samse ${hg38_bwa_index}/genome.fa - $raw_data_path/${base}.fq.gz | \
    samtools view -Shb -@ ${threads} - | \
    samtools sort -@ ${threads} -m 1G -o ${base}.bam -
done

exit(0)
#---
# merge two replicates in ChIP-seq data
# SRR1029001.bam & SRR1029002.bam
# are pseudo file names for next convenience
# and avoid to erase the alignment results.
# SRR1029001.bam will be heart_11.5.treat
# SRR1029002.bam will be heart_11.5.control
#---

##samtools merge -@ ${threads} -h SRR1029874.bam -O BAM SRR1029001.bam SRR1029874.bam SRR1029876.bam
##samtools merge -@ ${threads} -h SRR1029875.bam -O BAM SRR1029002.bam SRR1029875.bam SRR1029877.bam
##
##EOF
##
##echo 'fish step chunck'
##
##
##heart_whole_bams=($(find . -maxdepth 1 -type f -name "*.bam"))
##
##index_array=($(seq 22 1 33))
##index_array+=(48 49)
##index_num=${#index_array[@]}
##
###---
### original protocol and parameter setting
### is from the published paper
### here is the param excerpted from the 
### paper: S1, pdf = p12
### MACS call:
### macs14 -t chip:bam --control = input:bam --name = chip_output \
### --format =BAM --gsize =mm --tsize = 50 --bw = 300 --mfold \
### = 10; 30 --nolambda --nomodel --shiftsize = 150 -p 0:00001
###---
##
##for (( i=0; i<$((index_num)); i++ ));
##do
##    index=${index_array[$i]}
##    treat=${heart_whole_bams[$index]}
##    i=$((i+1))
##    index=${index_array[$i]}
##    control=${heart_whole_bams[$index]}
##    base=${treat%.bam}
##    macs2 callpeak --treatment ${treat} --control ${control} \
##          --format BAM --gsize mm --name ${base} --bdg --qvalue 0.01
##done
##
##source deactivate macs2
##