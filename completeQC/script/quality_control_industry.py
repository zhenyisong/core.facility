#!/usr/bin/env python

#---
# @author Yisong Zhen
# @since  2018-01-24
# @update 2018-05-24
#---

#---
# project aim:
# quality control pipeline in python3 version
# this package or script is designed for usage in 
# the core facility at GuoZhong.
#---

"""
the script setting path: 
if I used the cProfile module (with -m cProfile -s cumulative), 
the error will be thrown out.

Python multiprocessing pickling error
see here:
https://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error

python -m cProfile -s cumulative /wa/zhenyisong/sourcecode/core.facility/qualiControl/quality_control_industry.py -n '.downsample.fq.gz'
python quality_control_industry.py -n '.downsample.fq.gz' -g 'hg38' -l 'PE' -s 'RF'
python  /wa/zhenyisong/sourcecode/core.facility/qualiControl/code/quality_control_industry.py \
-n 'fastq.gz' -g 'mm10' -l 'PE' -s 'NONE' -w /home/zhenyisong/data/temp \
-d /home/zhenyisong/data/cardiodata/test/SRP109298


python  /wa/zhenyisong/sourcecode/core.facility/qualiControl/code/quality_control_industry.py \
   -g 'mm10' -l 'SE' -s 'NONE'  -u True

clean the all the output

rm -rf *.bam *.bai *.picard *.pdf fastqc* *.txt
"""


#---
# test set collections
# PE + stranded:     GSE85727     PMID: 27591185 SRP08239  Mouse
#                    GSE106688                   SRP124631 Human
# PE + Non-stranded: GSE10088     PMID: 28724790 SRP109298 Mouse
#                    GSE64403     PMID: 25477501 SRP051406 Mouse
# SE + stranded:     GSE77795     PMID: 27716916 SRP069873 Mouse
# SE + Non-stranded: GSE67889     PMID: 26641715 SRP057171
#
# script
# working snnipet without using the shell function
# I will optimize the following code someday.
# 
# Old way:
# to test this script is good,
# I downsample the data to 1% of the raw data
# source activate biotools
# cp /home/zhenyisong/data/results/chenlab/xiaoning/data/A_1/*.fq.gz /home/zhenyisong/data/temp/test
# seqtk sample -s 100 A_1_R1.fq.gz 0.01 | gzip - > A_1_R1.downsample.fq.gz
# seqtk sample -s 100 A_1_R2.fq.gz 0.01 | gzip - > A_1_R2.downsample.fq.gz
# this  will speed up the developemnt of the pyton3 QC pipeline
#
# https://www.biostars.org/p/121336/
# BBMap package? I did not try this method.
# Question: Select sequences from fastq.gz file
# several downsample methods to perform the same function above.
#---

"""
source activate biotools

#PE
 find ./SRP074376 -name '*.sra' | \
 xargs -P 3 -n 1 -I{} fastq-dump --outdir test/SRP074376 \
 --gzip --split-files --skip-technical {}
#SE
 find ./SRP069873  -name '*.sra' | \
 xargs -P 3 -n 1 -I{} fastq-dump --outdir test/SRP069873  \
 --gzip --skip-technical {}
 
find -name '*.fastq.gz' | xargs -I{} basename {} '.fastq.gz' | \
xargs -P 4 -n 1 -I{} sh -c 'seqtk sample -s100 {}.fastq.gz 0.01 | \
gzip -f - > {}.downsample.fq.gz'




#conda install parallel
# I tried, but failed due to the path transfer was incorrect.
#find  -name '*.fastq.gz' | \
#parallel --max-proc=3 'seqtk sample -s100 {/.} 300 | gzip - > {/.}.downsample.fq.gz'

"""
import sys
import completeQC.qc_modules as QC
import argparse
import timeit
import os




#---
# python style PEP8
# style guide for Python Code:
# constant definition use the UPPER CASE
#---

# claim default setting
# and initilize the variables needed in 
# the current script. These variables are 
# required by the PICARD modules. 
#---

STRANDNESS         = None
BWA_INDEX_PATH     = None
REFERENCE_GENOME   = None
THREADS            = 3

#--- default CONSTANT end


#---
# how to get the outside parameter from script 
# command line
# simple argparse example wanted: 1 argument, 3 results
# https://stackoverflow.com/questions/7427101/simple-argparse-example-wanted-1-argument-3-results
# in the following setting
# I created 6 parameters which can be recieved from the srcipt inputs
#
# Writing a help for python script
# https://stackoverflow.com/questions/9037828/writing-a-help-for-python-script
#---


param_parser = argparse.ArgumentParser(
                  prog        = ''' QC hunter''',
                  usage       = ''' python quality_control_industry.py -n '.downsample.fq.gz' 
                                    -g 'hg38' -l 'PE' -s 'RF' ''',
                  description = ''' This is a NGS Quality Control (QC) pipeline
                                    for the internal and academic usage only. Commercial
                                    permission should be solicitued from
                                    the State Laboratory for Cardiovascular Diseases, China. 
                                     ''',
                  epilog      = '''The script author is @zhenyisong;
                                   zhenyisong@fuwaihospital.org. Office-B201
                                   Tel: 86-01-60866301''')
param_parser.add_argument( '-q', '--QC-type', default = 'mRNA', 
                           choices = [ 'mRNA', 'miRNA', 'lncRNA',
                                       'ChIPseq','DNAseq'],
                           help    = """ this will choose the QC type;
                                         different NGS data will use different QC
                                         strategy. """)

param_parser.add_argument( '-s', '--strandness', default = 'NONE',
                           choices = [ 'NONE', 'FR','RF'],
                           help = """ this will set the strandness 
                                      in mNRA or lncRNA library
                                      construction method, which 
                                      is valuable to mapping and QC strategy.
                                      If the library uses the first strand to construct
                                      cDNA library, user should use the ___, and if using
                                      the second strand to construct library, user should 
                                      use the _____ """ )              
param_parser.add_argument( '-d', '--data-path', 
                           default   = os.getcwd(), 
                           required  = False,
                           help      = """ this will set the raw data path, 
                                           which the script can
                                           read sequencing data from this path""" )
param_parser.add_argument( '-g', '--genome-build', 
                           default = 'mm10', 
                           choices = [ 'hg38', 'hg19', 'mm10',
                                       'mm9','rn6'],
                           help    = """ this will set the genome version to quality, 
                                         control the data. For example, if the raw data
                                         are from human sample, then the most recent 
                                         version hg38 is recommended. If they are 
                                         from mouse or rat, the mm10
                                         and rn6 is recommended respectively """ )
param_parser.add_argument( '-l', '--library-model', 
                           default = 'PE', 
                           choices = [ 'PE','SE','MA'],
                           help    = """ this will set the read whether is paired or not, 
                                         single end - SE, paired end - PE. This parameter
                                         is derived from raw data sequencing model """ )
param_parser.add_argument( '-n', '--name-pattern', 
                           default   = '.fq.gz',
                           help      = """ this will set the read files suffix. The 
                                           script will find the file with secified
                                           suffix and use them as the raw data input """ )
param_parser.add_argument( '-w', '--working-path', 
                           default  = './', 
                           required = False,
                           help     = """ this will set the output directory and the Python3
                                          script working directory. This dir will have 
                                          writable and readable permission """ )
param_parser.add_argument( '-t', '--threads', 
                           default   = 3, 
                           type      = int, 
                           required  = False,
                           help      = """ this will set the thread number which is used in 
                                           fastqc module and bwa module """ )
param_parser.add_argument( '-a', '--aligner', 
                           default = 'BWA', 
                           choices = ['BWA','HISAT2'],
                           help    = """ this will set the alignment algorithm used in 
                                         the alignment procedure when to genenrate BAM files """ )
param_parser.add_argument( '-u', '--debug', 
                           default = False, 
                           type    = bool,
                           choices = [True,False],
                           help    = """ this will use the sample set to simulate the  
                                         QC pipeline results """ )
   
start_time   = timeit.default_timer()
whole_params = vars(param_parser.parse_args())

QC.perform_mRNA_QCtask(whole_params)

stop_time    = timeit.default_timer()

QC.compute_running_time(start_time, stop_time)

sys.exit(0)

