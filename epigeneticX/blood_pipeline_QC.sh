#!/bin/bash

#---
#@author Yisong Zhen
#@since  2018-03-26
#update  2018-04-02
#---

# get the Mol. Cell data by using aspera
# GSE31332
# 
# tar xvzf igenomes/Drosophila_melanogaster_UCSC_dm6.tar.gz
# 
# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/blood_pipeline_QC.sh
# nohup bash /wa/zhenyisong/sourcecode/core.facility/epigeneticX/blood_pipeline_QC.sh &


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


#---
#find /home/zhenyisong/data/cardiodata/SRP007863  -name '*.sra' | \
#xargs -P 3 -n 1 -I{} fastq-dump \
#--outdir /home/zhenyisong/data/results/chenlab/songli/pipelineQC \
#--gzip --skip-technical {}
#
#  gunzip SRR360699.fastq.gz  SRR360700.fastq.gz SRR360701.fastq.gz
#---

QC_pipeline_raw_data='/home/zhenyisong/data/results/chenlab/songli/pipelineQC/rawdata'
mapping_bwa_results='/home/zhenyisong/data/results/chenlab/songli/pipelineQC/bwa'
fly10_bwa_index='/wa/zhenyisong/reference/Drosophila_melanogaster/UCSC/dm6/Sequence/BWAIndex'
fly10_dict_index='/wa/zhenyisong/reference/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa'
fly10_bowtie_1_index='/home/zhenyisong/data/reference/Drosophila_melanogaster/UCSC/dm6/Sequence/BowtieIndex'
read_length=36

#---
# wget -np -nd -r http://hgdownload.cse.ucsc.edu/goldenPath/dm6/bigZips/dm6.chrom.sizes
# or use the command provided by UCSC
# 
#---
dm6_chrom_sizes='/home/zhenyisong/data/results/chenlab/songli/chirpseq/dm6.chrom.sizes'
sam2bedGraph='/home/zhenyisong/data/results/chenlab/songli/chirpseq/sam2bedGraph_norm.py'
merge2bedGraph='/home/zhenyisong/data/results/chenlab/songli/chirpseq/combine_two_bedGraph.pl'
bedGraph2sam='/home/zhenyisong/data/results/chenlab/songli/chirpseq/bedGraph2sam.pl'
peakCorrelation='/home/zhenyisong/data/results/chenlab/songli/chirpseq/peak_correlation.pl'


if [ -d "${mapping_bwa_results}" ];then
    rm -rf ${mapping_bwa_results}
    mkdir -p ${mapping_bwa_results}
    cd ${mapping_bwa_results}
fi

if [ ! -f 'genome.*.ebwt' ];then
    ln -s "${fly10_bowtie_1_index}"/genome.* ./
fi

# bowtie 1 only accept decompressed files
#---
fly_roX2_fastq=(SRR360699.fastq SRR360700.fastq SRR360701.fastq)

file_number=${#fly_roX2_fastq[@]}

for (( i=0; i<$((file_number)); i++ ));
do
    filename=${fly_roX2_fastq[$i]}
    base=`basename ${filename}`
    base=${base%.fastq}
    bowtie -m 1 -p 4 -S -q genome ${QC_pipeline_raw_data}/${filename} > ${base}.sam
done



python ${sam2bedGraph} SRR360699.sam even.molCell.bedGraph
python ${sam2bedGraph} SRR360700.sam odd.molCell.bedGraph
python ${sam2bedGraph} SRR360701.sam control.molCell.bedGraph

perl ${merge2bedGraph} ${read_length} ${dm6_chrom_sizes} \
                                      even.molCell.bedGraph \
                                      odd.molCell.bedGraph \
                                      merge.molCell.wig

wigToBigWig merge.molCell.wig ${dm6_chrom_sizes} merge.molCell.bw
bigWigToBedGraph merge.molCell.bw merge.molCell.bedGraph

perl ${bedGraph2sam} ${dm6_chrom_sizes} merge.molCell.bedGraph merge.molCell.sam

bedGraphToBigWig even.molCell.bedGraph ${dm6_chrom_sizes} even.molCell.bw
bedGraphToBigWig odd.molCell.bedGraph  ${dm6_chrom_sizes} odd.molCell.bw
bedGraphToBigWig control.molCell.bedGraph ${dm6_chrom_sizes} control.molCell.bw

source deactivate macs2

macs14 -t merge.molCell.sam -c SRR360701.sam -f SAM -n merge --bw 200 -m 10,50
perl ${peakCorrelation} merge_peaks.xls \
                        even.molCell.bedGraph \
                        odd.molCell.bedGraph \
                        merge.molCell.bedGraph