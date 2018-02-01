#---
# project aim:
# quality control pipeline in python3 version
# this package or script is designed for usage in 
# the core facility at GuoZhong.
#---

#---
# @author Yisong Zhen
# @since  2018-01-24
# @update 2018-02-01
#---

import os
import re
import glob
import tempfile
from plumbum import local, FG, BG
from plumbum.cmd import cut, rm

#---
# python style PEP8
# style guide for Python Code:
# constant definition use the UPPER CASE
#---

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

HISAT2_INDEX_MM10_PATH   = '/wa/zhenyisong/reference/index/mm10'
BWA_INDEX_MM10_PATH      = '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa'


#---
# define the program threads
# this parameter is required by aligner program BWA or HISAT2
#---
THREADS = 6


#---
# now the script get the input parameters
# which are specified by the end user.
#---

raw_data_path    = '/home/zhenyisong/data/results/chenlab/xiaoning/data'
raw_data_pattern = '/home/zhenyisong/data/temp/test'
working_dir      = '/home/zhenyisong/data/temp/test'
os.chdir(working_dir)


#---
# Linux envs uses conda to manage the biological
# softwares, which mush pre-install these module
# to perform the task implemented in this python
# script.
# aim -- use the Plumbum python package to mimic
#        the linux shell command,
#     -- define these linux commands in shell way
#---

hisat2   = local['hisat2']
samtools = local['samtools']
gatk4    = local['gatk-launch']
picard   = local['picard']
fastqc   = local['fastqc']
bwa      = local['bwa']


'''
read_1_files      = []
read_2_files      = []
ending_pattern_1  = '_R1.downsample.fq.gz'
ending_pattern_2  = '_R2.downsample.fq.gz'
for file in glob.glob(raw_data_pattern, recursive = True):
	if file.endswith(ending_pattern_1):
		read_1_files.append(file)
	elif file.endswith(ending_pattern_2):
		read_2_files.append(file)
'''

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

'''
@parameters needed
   1. threads: for paraelle computition
   2. reads file names; two model PE or SE model
   3. reference geneome location
   4. out_put dir
@function
   to perform the BWA alignment and output the
   coordination sorted BAM format alignment files

@test
   run_BWA_aligner(read_1_files[0], read_2_files[0]) 
@return
   the read1 and read2 file names (absolute paths)
   in tuple.

'''

def run_BWA_aligner( read1, read2,
                     bwa_index_file  = BWA_INDEX_MM10_PATH,
                     sam_index_file  = MM10_UCSC_GENOME,
                     threads     = THREADS ):
    basename = get_basename(read1, '_R1.downsample.fq.gz')
    run_bwa  = ( 
         bwa[ 
             'mem',
             '-M',
             '-t',threads,
             bwa_index_file,
             read1, read2
         ] | samtools[
             'view',
             '-Sbh',
             '-@',threads,
             '-O', 'BAM',
             '-T', sam_index_file
         ] | picard[
             'SortSam', 
             'INPUT=','/dev/stdin',
             'OUTPUT=', basename + '.bam',
             'SORT_ORDER=','coordinate'
         ]
      )
    run_bwa()
    return (read1,read2)

'''
@parameters needed
   1. threads: for paraelle computition
   2. reads file names; two model PE or SE model
   3. reference geneome location
   4. out_put dir
@function
   to perform the HISAT2 alignment and output the
   coordination sorted BAM format alignment files

@test
   run_HISAT2_aligner(read_1_files[0], read_2_files[0]) 
@return
   the read1 and read2 file names (absolute paths)
   in tuple.

'''
def run_HISAT2_aligner( read1, read2,
                        hisat2_index_file = HISAT2_INDEX_MM10_PATH,
                        sam_index_file    = MM10_UCSC_GENOME,
                        threads           = threads ):
    hisat2_run = (
          hisat2[
              '-p', threads,
              '--dta',
              '--fr', 
              '-x', hisat2_index_file, 
              '-1', read1, '-2', read2,
          ] | samtools[
              'view',
              '-Sbh',
              '-@',threads,
              '-O', 'BAM',
              '-T', sam_index_file
          ] | picard[
              'SortSam', 
              'INPUT=','/dev/stdin',
              'OUTPUT=', output_filename,
              'SORT_ORDER=','coordinate'
          ]
        )

    hisat2_run()
    return(read1, read2)


'''
@aim 
    get the bam index used the GATK4 program
    this will speed the program scanning of the 
    bam result.
@parameters
    input_file (Stirng):    the bam file which need to be indexed
    output_file(String):   output of the indexed bam file.
@return
    bam indexed file name. String

'''

def build_BAM_index(input_file = , output_file = ):
    
    run_build_index = ( gatk4[ 
                          'BuildBamIndex', 
                          '--INPUT',  input_file,
                          '--OUTPUT', output_file] )
    run_build_index()
    return output_file

'''
@aim 
    get the base name of the sequncing raw data for
    interna usage as the sample names.
@parameters
    fullname:         file full name
    ending_pattern:   the read data ending pattern which 
                      discriminate the raw data from other 
                      non-related files
@return
    the file base-name which may be used as the sample name
    to trace the QC results with the corresponding raw data

'''

def get_basename(fullname, ending_pattern = 'fq.gz' ):
    basename = os.path.basename(fullname)
    basename = re.sub(ending_pattern,'', basename)
    return basename

'''
@aim 


@parameters
    sample_name:      get the sample name (with base + '.bam') sample name
                      was setted to be the basename free of ending pattern.
    ref_genome:       the location of reference genome of selected animal 
                      model. this will be the
    ref_flat:         the ref_flat file generated according to the suggestion
                      by BioStar post.
    strandness:       the input parameter transfered from the python script

@return
    the QC result by PICARD modules.

'''
def run_PICARD_QCmodules( sample_name = ,
                          ref_genome = MM10_UCSC_GENOME,
                          ref_flat   = REFFLAT_MM10_UCSC_PICARD,
                          strandness = STRANDNESS):

    picard_run = ( picard[ 'CollectRnaSeqMetrics',
                           'REF_FLAT=', ref_flat,
                           'RIBOSOMAL_INTERVALS=',TEMP_FILE_NAME,
                           'STRAND_SPECIFICITY=',strandness,
                           'CHART_OUTPUT=','null',
                           'METRIC_ACCUMULATION_LEVEL=', 'ALL_READS',
                           'INPUT=',  base + '.bam',
                           'OUTPUT=', base + '.CollectRnaSeqMetrics.picard',
                           'ASSUME_SORTED=','true'] )
    picard_run()

    picard_run = ( picard[ 'CollectAlignmentSummaryMetrics',
                           'REFERENCE_SEQUENCE=', ref_genome,
                           'INPUT=',  base + '.bam',
                           'OUTPUT=', base + '.CollectAlignmentSummaryMetrics.picard',
                           'EXPECTED_PAIR_ORIENTATIONS=','null']
                  )
    picard_run()

    picard_run = ( picard[ 'CollectInsertSizeMetrics',
                           'INPUT=',  base + '.bam',
                           'OUTPUT=', base + '.bam',
                           'HISTOGRAM_FILE=','tmp.foool.pdf',
                           'MINIMUM_PCT=', 0.05] )
    picard_run()
    
    
    remove_file = (rm[TEMP_FILE_NAME])
    remove_file()


'''
@aims
    the function parse the data path and fetch all
    raw data required to perform QC checking
    raw read data are using the relative or absolute
    path.
@parameters
    raw_data_path: the raw illumina reads sequencing results
    end_pattern:   the read data ending pattern which 
                   discriminate them from other non-related
                   files
@return
    the function return the file lists

'''
def get_READSeq_files(raw_data_path, ending_pattern = '.fq.gz'):
    raw_data_path = raw_data_path + '/**/*' + ending_pattern
    files         = []
    for file in glob.glob(raw_data_path, recursive = True):
        files.append(file)
    return files

'''
@aim:
    to get the Pair-end sequencing file list
    read1 and read2 file list.
@parameters
    files: the raw illumina reads file list

@return
    the splitted files in READ1 and READ2 format
    list

'''
def split_PairEnd_files(files):
    read1_list = files[1::2]
    read2_list = files[::2]
    return read1_list, read2_list


def run_FASTQC( file, threads = THREADS,
                output_dir = 'fastqc.results'):
    os.makedirs(fastqc_dir, exist_ok = True)
    run_fastqc = ( fastqc[ '-o', fastqc_dir,
                           '-f', 'fastq',
                           '-t', threads,
                           file] )
    run_fastqc()

def get_RIBO_file( file, ribo_annotation = RIBO_INTERVAL_LIST_MM10_PICARD):
    TEMP_FILE      = tempfile.NamedTemporaryFile(dir = os.getcwd())
    TEMP_FILE_NAME = TEMP_FILE.name
    TEMP_FILE.close()
    file_basename  =
    get_samtool_header = ( samtools[ 'view',
                                     '-H', base + '.bam',
                                     '-o', TEMP_FILE_NAME] )
    get_samtool_header()
    
    get_ribo_file  = ( cut[ '-s',
                            '-f', '1,4,5,7,9',
                            ribo_annotation] >> TEMP_FILE_NAME )
    get_ribo_file()

def choose_ANNOTATION():


