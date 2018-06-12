#!/bin/bash
# @author Yisong Zhen
# @since  2018-06-11
# @update 2018-06-11
#---



WGBS_index='/wa/zhenyisong/reference/WGBS/mouse'
#---
# conda install bismark
# cd ${WGBS_index}
# bismark_genome_preparation --bowtie2 ./
#---
bismark --multicore 3 --bowtie2 --bam ${WGBS_index}