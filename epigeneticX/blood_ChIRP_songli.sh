#!/bin/bash
# @author  Yisong Zhen
# @since   2018-03-22
# @update  2018-03-28
#---


# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/blood_ChIRP_songli.sh
# nohup bash /wa/zhenyisong/sourcecode/core.facility/epigeneticX/blood_ChIRP_songli.sh &


#----
# HPC parameters for Sun Grid
#$ -N Yisong.MACS2
#$ -S /bin/bash
#$ -w e
#$ -wd /home/zhenyisong/data/results/chenlab/songli
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

#---
# 
# the original data were uploaded from songli
# she has a project on CHiRP.
# I have already checked the data fingerprints
# using the command
# md5sum -c CleanFq_md5sum.TXT
# and all data are intact and in good state
#--


threads=4


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
hg38_dict_index='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'


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

all_raw_data_read_1=($raw_data_path/*_1.fq.gz)
all_raw_data_read_2=($raw_data_path/*_2.fq.gz)
file_number=${#all_raw_data_read_1[@]}


echo 'my own protocol step for songli data'
: << 'EOF'



for (( i=0; i<$((file_number)); i++ ));
do
    filename=${all_raw_data_read_1[$i]}
    base=`basename ${filename}`
    base=${base%_1.fq.gz}
    #---
    # I used the -M parameter to cater for the PICARD pipeline
    # see: https://www.biostars.org/p/234768/
    # just keep the default setting
    #---
    bwa mem -M -t ${threads} ${hg38_bwa_index}/genome.fa \
    $raw_data_path/${base}_1.fq.gz $raw_data_path/${base}_2.fq.gz| \
    picard SortSam  INPUT=/dev/stdin OUTPUT="${base}.bam" SORT_ORDER=coordinate
done

fastqc -f fastq -t ${threads} -q -o ${fastqc_results} $raw_data_path/*.fq.gz


EOF

for (( i=0; i<$((file_number)); i++ ));
do
    filename=${all_raw_data_read_1[$i]}
    base=`basename ${filename}`
    base=${base%_1.fq.gz}
    picard CollectAlignmentSummaryMetrics R="${hg38_dict_index}" \
                                          I="${base}.bam" \
                                          O="${base}.txt"
    picard CollectInsertSizeMetrics I="${base}.bam" \
           O="${base}.insert_size_metrics.txt" \
           H="${base}.insert_size_histogram.pdf" M=0.5
done

cd ${mapping_bwa_results}

treat_bam_files=(Even_clean.bam ODD_clean.bam)
control_bam_files=(Input_clean.bam Input_clean.bam)
index_num=${#control_bam_files[@]}


for (( i=0; i<$((index_num)); i++ ));
do
    treat=${treat_bam_files[$i]}
    control=${control_bam_files[$i]}
    base=${treat%.bam}
    macs2 callpeak --treatment ${treat} --control ${control} \
          --format BAMPE --gsize hs --name ${base} --bdg --qvalue 0.01
done

samtools merge -@ 4 -f blood.full.bam Even_clean.bam ODD_clean.bam
macs2 callpeak --treatment blood.full.bam --control Input_clean.bam \
          --format BAMPE --gsize hs --name blood --bdg --qvalue 0.01


# extra replication study to compare with t
# the nova analysis result
# I did not replicate the result
# by S-roX2 ChIRP data and have to add extra
# QC steps.
#--
nova_path='/home/zhenyisong/data/results/chenlab/songli/BAM'
nova_bams=($nova_path/*.bam)
nova_file_number=${#nova_bams[@]}

for (( i=0; i<$((nova_file_number)); i++ ));
do
    filename=${nova_bams[$i]}
    base=`basename ${filename}`
    base=${base%.bam}
    picard CollectAlignmentSummaryMetrics R="${hg38_dict_index}" \
                                          I="${base}.bam" \
                                          O="${base}.nova.txt"
    picard CollectInsertSizeMetrics I="${base}.bam" \
           O="${base}.insert_size_metrics.nova.txt" \
           H="${base}.insert_size_histogram.nova.pdf" M=0.5
done


treat_bam_files=(Even.bam ODD.bam)
control_bam_files=(Input.bam Input.bam)
index_num=${#control_bam_files[@]}
for (( i=0; i<$((nova_file_number)); i++ ));
do
    treat=${treat_bam_files[$i]}
    control=${control_bam_files[$i]}
    base=${treat%.bam}
    macs2 callpeak --treatment ${nova_path}/${treat} --control ${nova_path}/${control} \
          --format BAMPE --gsize hs --name "${base}.nova" --bdg --qvalue 0.01
done

#---
# macs2 only accepts sorted bam files. Please see the github feedback by Tao Liu.
#---
samtools merge -@ 4 -f blood.nova_full.bam ${nova_path}/Even.bam ${nova_path}/ODD.bam
picard SortSam  INPUT='blood.nova_full.bam' OUTPUT="blood.nova_sort.bam" SORT_ORDER=coordinate
mv blood.nova_full.bam blood.nova_merge.bam
picard BuildBamIndex INPUT='blood.nova_merge.bam'
macs2 callpeak --treatment blood.nova_merge.bam --control ${nova_path}/Input.bam \
               --format BAMPE --gsize hs --name blood_nova --bdg --qvalue 0.01

bamCoverage --bam ${nova_path}/Even.bam --outFileFormat bedgraph --outFileName even.nova.bedgraph
bamCoverage --bam ${nova_path}/ODD.bam  --outFileFormat bedgraph --outFileName odd.nova.bedgraph
bamCoverage --bam blood.nova_merge.bam  --outFileFormat bedgraph --outFileName merge.nova.bedgraph

#---
# peak_correlation.pl
# this script is from this Github account;
# https://github.com/bdo311
# and I mirorred his codes on ChIRP
# https://github.com/bdo311/chirpseq-analysis
# git clone git@github.com:bdo311/chirpseq-analysis.git
#---
chirp_correlation='/home/zhenyisong/data/results/chenlab/songli/chirpseq-analysis/peak_correlation.pl'
ln -s ${chirp_correlation} ./
perl peak_correlation.pl blood_nova_peaks.xls even.nova.bedgraph odd.nova.bedgraph merge.nova.bedgraph

source deactivate macs2
