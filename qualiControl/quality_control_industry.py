#---
# project aim:
# quality control pipeline in python3 version
# this package or script is designed for usage in 
# core facility
#---

#---
# @author Yisong Zhen
# @since  2018-01-24
# @update 2018-01-30
#---

import os
import re
import glob
import tempfile
from plumbum import local, FG, BG
from plumbum.cmd import cut, rm

# python style:
# style guide for Python Code:
# constant definition use the UPPER CASE
#

#---
# annotation and genome files location
# location in linux platform is fixed
# the reference genome sets for various spieces
# were downloaded from iGenomes
#---
MM10_UCSC_GENOME = '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
HG38_UCSC_GENOME = '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
MM10_UCSC_GTF    = '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
HG38_UCSC_GTF    = '/home/zhenyisong/data/reference/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'

#---
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

RRNA_HG38_RSEQC                = '/wa/zhenyisong/reference/annotation/RSeQC/hg38_rRNA.bed'
HOUSE_KEEPING_GENES_HG38_RSEQC = '/wa/zhenyisong/reference/annotation/RSeQC/hg38.HouseKeepingGenes.bed'
BASIC_GENES_GENCODE_HG38_RSEQC = '/wa/zhenyisong/reference/annotation/RSeQC/hg38_GENCODE_v24_basic.bed'

RRNA_MM10_RSEQC                = '/wa/zhenyisong/reference/annotation/RSeQC/mm10_rRNA.bed'
HOUSE_KEEPING_GENES_MM10_RSEQC = '/wa/zhenyisong/reference/annotation/RSeQC/mm10.HouseKeepingGenes.bed'
BASIC_GENES_GENCODE_MM10_RSEQC = '/wa/zhenyisong/reference/annotation/RSeQC/mm10_GENCODE_VM11_basic.bed'


#---
# picard annotation files
# see more help on picard
# http://broadinstitute.github.io/picard/command-line-overview.html
#
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
#
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

REFFLAT_MM10_UCSC_PICARD       = '/wa/zhenyisong/reference/annotation/picard/refFlat_mm10.txt'
RIBO_INTERVAL_LIST_MM10_PICARD = '/wa/zhenyisong/reference/annotation/picard/mm10_ribosome_interval_list.txt'

#---
# hisat2 index file location
# and other required annotation files
#---

HISAT2_INDEX_PATH = '/home/zhenyisong/data/reference/index/mm10'


# define the program threads
threads = 6


# now the script get the input parameters
#

raw_data_path    = '/home/zhenyisong/data/results/chenlab/xiaoning/data'
raw_data_pattern = '/home/zhenyisong/data/temp/test/**/*.fq.gz'
working_dir      = '/home/zhenyisong/data/temp/test'
os.chdir(working_dir)


#---
# module I
# aim -- generate bam file from raw data,
#     -- genome read mapping
#---

hisat2   = local['hisat2']
samtools = local['samtools']
gatk4    = local['gatk-launch']
picard   = local['picard']
fastqc   = local['fastqc']



## if unpiared the data, -U parameter will be used
##shopt -s nullglob

read_1_files      = []
read_2_files      = []
ending_pattern_1  = '_R1.downsample.fq.gz'
ending_pattern_2  = '_R2.downsample.fq.gz'
for file in glob.glob(raw_data_pattern, recursive = True):
	if file.endswith(ending_pattern_1):
		read_1_files.append(file)
	elif file.endswith(ending_pattern_2):
		read_2_files.append(file)

#---
# to test this script is good,
# I downsample the data to 2% of the raw data
# cp /home/zhenyisong/data/results/chenlab/xiaoning/data/A_1/*.fq.gz /home/zhenyisong/data/temp/test
# seqtk sample -s 100 A_1_R1.fq.gz 0.2 | gzip - > A_1_R1.downsample.fq.gz
# seqtk sample -s 100 A_1_R2.fq.gz 0.2 | gzip - > A_1_R2.downsample.fq.gz
# this  will speed up the developemnt of the pyton3 QC pipeline
#---


for i in range(len(read_1_files)):
    R1   = read_1_files[i]
    R2   = read_2_files[i]
    base = os.path.basename(R1)
    base = re.sub(ending_pattern_1,'', base)
    output_filename    = base + '.bam'
    output_indexname   = base + '.bam.bai'
    pipeline = (
          hisat2[
              '-p', threads,
              '--dta',
              '--fr', 
              '-x', HISAT2_INDEX_PATH, 
              '-1', R1, '-2', R2,
          ] | samtools[
              'view',
              '-Sbh',
              '-@',threads,
              '-O', 'BAM',
              '-T', MM10_UCSC_GENOME
          ] | picard[
              'SortSam', 
              'INPUT=','/dev/stdin',
              'OUTPUT=', output_filename,
              'SORT_ORDER=','coordinate'
          ]
        )
    pipeline()
    build_index = ( gatk4[ 'BuildBamIndex', 
                         '--INPUT',output_filename,
                         '--OUTPUT', output_indexname] )
    build_index()

fastqc_dir = 'fastqc.results'
os.makedirs(fastqc_dir, exist_ok = True)

for i in range(len(read_1_files)):
    R1   = read_1_files[i]
    R2   = read_2_files[i]
    fastqc_run = ( fastqc[ '-o', fastqc_dir,
                           '-f', 'fastq',
                           '-t', threads,
                           R1, R2 ] )
    fastqc_run()


#---
# module 3
# aim -- using the picard procedure to infer the QC
#     -- rRNA percentage etc.
#---
STRANDNESS = 'NONE'



for i in range(len(read_1_files)):
    R1   = read_1_files[i]
    R2   = read_2_files[i]
    base = os.path.basename(R1)
    base = re.sub(ending_pattern_1,'', base)
    TEMP_FILE      = tempfile.NamedTemporaryFile(dir = os.getcwd())
    TEMP_FILE_NAME = TEMP_FILE.name
    TEMP_FILE.close()
    get_samtool_header = ( samtools[ 'view',
                                     '-H', base + '.bam',
                                     '-o', TEMP_FILE_NAME] )
    get_samtool_header()
    
    get_ribo_file  = ( cut[ '-s',
                            '-f', '1,4,5,7,9',
                            RIBO_INTERVAL_LIST_MM10_PICARD] >> TEMP_FILE_NAME )
    get_ribo_file()
    '''
    picard_run = ( picard[ 'CollectRnaSeqMetrics',
                           'REF_FLAT=', REFFLAT_MM10_UCSC_PICARD,
                           'RIBOSOMAL_INTERVALS=',TEMP_FILE_NAME,
                           'STRAND_SPECIFICITY=',STRANDNESS,
                           'CHART_OUTPUT=','null',
                           'METRIC_ACCUMULATION_LEVEL=', 'ALL_READS',
                           'INPUT=',  base + '.bam',
                           'OUTPUT=', base + '.CollectRnaSeqMetrics.picard',
                           'ASSUME_SORTED=','true'] )
    picard_run()

    picard_run = ( picard[ 'CollectAlignmentSummaryMetrics',
                           'REFERENCE_SEQUENCE=', MM10_UCSC_GENOME,
                           'INPUT=',  base + '.bam',
                           'OUTPUT=', base + '.CollectAlignmentSummaryMetrics.picard',
                           'EXPECTED_PAIR_ORIENTATIONS=','null']
                  )
    picard_run()
    '''
    picard_run = ( picard[ 'CollectInsertSizeMetrics',
                           'INPUT=',  base + '.bam',
                           'OUTPUT=', base + '.bam',
                           'HISTOGRAM_FILE=','tmp.foool.pdf',
                           'MINIMUM_PCT=', 0.05] )
    picard_run()
    
    
    remove_file = (rm[TEMP_FILE_NAME])
    remove_file()


#source deactivate macs2

#multiqc ./