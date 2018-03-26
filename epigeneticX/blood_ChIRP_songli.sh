#!/bin/bash
# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-03-26
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
# and CONSTANT setting
#---

mapping_bwa_results='/wa/zhenyisong/results/chenlab/songli/bwa'
fastqc_results='/wa/zhenyisong/results/chenlab/songli/bwa/fastqc'
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

if [ ! -d "${fastqc_results}" ]; then
    mkdir -p "${fastqc_results}"
else
    rm -rf "${fastqc_results}"
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
    #---
    # I used the -M parameter to cater for the PICARD pipeline
    # see: https://www.biostars.org/p/234768/
    # just keep the default setting
    #---
    bwa mem -M -t ${threads} ${hg38_bwa_index}/genome.fa $raw_data_path/${base}.fq.gz | \
    picard SortSam  INPUT=/dev/stdin OUTPUT={base}.bam SORT_ORDER=coordinate
done

fastqc -f fastq -t ${threads} -q -o ${fastqc_results} $raw_data_path/*.fq.gz

cd ${mapping_bwa_results}

treat_bam_files=(Even_clean_1.bam ODD_clean_1.bam Even_clean_2.bam ODD_clean_2.bam)
control_bam_files=(Input_clean_1.bam Input_clean_1.bam Input_clean_2.bam Input_clean_2.bam)
index_num=${#control_bam_files[@]}


for (( i=0; i<$((index_num)); i++ ));
do
    treat=${treat_bam_files[$i]}
    control=${control_bam_files[$i]}
    base=${treat%.bam}
    macs2 callpeak --treatment ${treat} --control ${control} \
          --format BAM --gsize hs --name ${base} --bdg --qvalue 0.01
done

source deactivate macs2
