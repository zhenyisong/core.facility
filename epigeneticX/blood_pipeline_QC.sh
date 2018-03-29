#!/bin/bash

#---
#@author Yisong Zhen
#@since  2018-03-26
#update  2018-03-28
#---

# get the Mol. Cell data by using aspera
# GSE31332
# 
# tar xvzf igenomes/Drosophila_melanogaster_UCSC_dm6.tar.gz
# 
# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/blood_pipeline_QC.sh


#----
# HPC parameters for Sun Grid
#$ -N songliQC
#$ -S /bin/bash
#$ -w e
#$ -wd /home/zhenyisong/data/results/chenlab/songli/pipelineQC
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
#---

#---
# I asked wenke, and confirmed the cluster management system
# I think we can search the net using google
# "sge qsub script example'
# follwoing script is runing on the macs2 env by conda management
#
#---

source ~/.bashrc
source ~/.bash_profile
source activate macs2
unset PYTHONPATH


#find /home/zhenyisong/data/cardiodata/SRP007863  -name '*.sra' | \
#xargs -P 3 -n 1 -I{} fastq-dump \
#--outdir /home/zhenyisong/data/results/chenlab/songli/pipelineQC \
#--gzip --skip-technical {}
#
#  gunzip SRR360699.fastq.gz  SRR360700.fastq.gz SRR360701.fastq.gz
#

songli_QCresults='/home/zhenyisong/data/results/chenlab/songli/pipelineQC'
mapping_bwa_results='/home/zhenyisong/data/results/chenlab/songli/pipelineQC/bwa'
fly10_bwa_index='/wa/zhenyisong/reference/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex'

cd "${songli_QCresults}"

if [ ! -f 'genome.fa.ann' ]; then
    ln -s "${fly10_bwa_index}"/genome.* ./
fi

fly_roX2_fastq=(SRR360699.fastq  SRR360700.fastq SRR360701.fastq)

file_number=${#fly_roX2_fastq[@]}

for (( i=0; i<$((file_number)); i++ ));
do
    filename=${fly_roX2_fastq[$i]}
    base=`basename ${filename}`
    base=${base%.fastq}
    bwa mem -M -t 4 genome.fa ${filename} | \
    picard SortSam INPUT=/dev/stdin OUTPUT=${base}.bam SORT_ORDER=coordinate
    picard BuildBamIndex INPUT=${base}.bam
done

treat_bam_files=(SRR360699.bam SRR360700.bam)
control_bam_files=(SRR360701.bam SRR360701.bam)
base_name=(even odd)
index_num=${#control_bam_files[@]}


for (( i=0; i<$((index_num)); i++ ));
do
    treat=${treat_bam_files[$i]}
    control=${control_bam_files[$i]}
    base=${base_name[$i]}
    macs2 callpeak --treatment ${treat} --control ${control} \
          --format BAM --gsize dm --name ${base} --bdg --qvalue 0.01
done


samtools merge -@ 4 -f full.bam SRR360700.bam SRR360699.bam
macs2 callpeak --treatment full.bam --control SRR360701.bam \
               --format BAM --gsize dm --name roX2 --bdg --qvalue 0.01

source deactivate macs2



