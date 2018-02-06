#---
# project aim:
# quality control pipeline in python3 version
# this package or script is designed for usage in 
# the core facility at GuoZhong.
#---


# the script setting path: 
#  python /wa/zhenyisong/sourcecode/core.facility/qualiControl/quality_control_industry.py
#---



#---
# @author Yisong Zhen
# @since  2018-01-24
# @update 2018-02-06
#---

import os
import sys
import re
import glob
import tempfile
import argparse
from plumbum import local, FG, BG
from plumbum.cmd import cut, rm

#---
# python style PEP8
# style guide for Python Code:
# constant definition use the UPPER CASE
#---


# how to get the outside parameter from script 
# command line
# simple argparse example wanted: 1 argument, 3 results
# https://stackoverflow.com/questions/7427101/simple-argparse-example-wanted-1-argument-3-results
#---

'''
param_parser = argparse.ArgumentParser()
param_parser.add_argument( '-f', '--file', default = 'mRNA', 
                           choices = [ 'mRNA', 'miRNA', 'lncRNA',
                                       'ChIPseq','DNAseq'])

param_parser.add_argument( '-s', '--strandness', default = "None",
                           help = """ this will set the strandness 
                                      in mNRA or lncRNA library
                                      construction method, which 
                                      is valuable to mapping and QC strategy """ )              
param_parser.add_argument( '-d', '--data-path', default = "None", required = True,
                           help = """ this will set the raw data path, 
                                      which user can
                                      read sequencing data from """ )
param_parser.add_argument( '-g', '--genome-version', default = 'mm10', required = True,
                           choices = [ 'hg38', 'hg19', 'mm10',
                                       'mm9','rn6'],
                           help = """ this will set the genome version to quality, 
                                      control the data """ )
param_parser.add_argument( '-l', '--library-model', default = 'PE', required = True,
                           choices = [ 'PE','SE','MA'],
                           help = """ this will set the read whether is paired or not, 
                                      single end - SE, paired end - PE """ )
param_dict = param_parser.parse_args()

'''

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

STRANDNESS = 'NONE'
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


#---
# to test this script is good,
# I downsample the data to 2% of the raw data
# cp /home/zhenyisong/data/results/chenlab/xiaoning/data/A_1/*.fq.gz /home/zhenyisong/data/temp/test
# seqtk sample -s 100 A_1_R1.fq.gz 0.2 | gzip - > A_1_R1.downsample.fq.gz
# seqtk sample -s 100 A_1_R2.fq.gz 0.2 | gzip - > A_1_R2.downsample.fq.gz
# this  will speed up the developemnt of the pyton3 QC pipeline
#---

data_PEfile_list         = get_READSeq_files( raw_data_pattern, 
                                              ending_pattern = '.downsample.fq.gz')
run_FASTQC(data_PEfile_list)
read1_list, read2_list   = split_PairEnd_files(data_PEfile_list )
run_BWA_aligner(data_PEfile_list[0], data_PEfile_list[1], ending_pattern = '.downsample.fq.gz')
run_HISAT2_aligner(data_PEfile_list[0], data_PEfile_list[1], ending_pattern = '.downsample.fq.gz')



build_BAM_index('A_1_R1.bam','A_1_R1.bam.bai')
run_PICARD_QC_modules('A_1_R1')


#source deactivate macs2

#multiqc ./

'''
@parameters needed
   1. threads: (integer) for paraelle computation. This param use 
               the default to THREADS previous defined.
   2. reads1 & read2: (String) the absolute read path in string format
                      file names; two model PE or SE model
   3. bwa_index_file (String) : the location of the indexed genome files.
   4. sam_index_file (String) : the required reference genome file used
                                by samtools index
   5. ending_pattern (String) : the raw data ending pattern, use of which
                                to extract sample (or base name)
@function
   to perform the BWA alignment and output the
   coordination sorted BAM format alignment files

@test
   run_BWA_aligner(read_1_files[0], read_2_files[0], ending_pattern = '.downsample.fq.gz') 
@return: 
   the read1 and read2 file names (absolute paths)
   in tuple.

'''

def run_BWA_aligner( read1, read2,
                     bwa_index_file  = BWA_INDEX_MM10_PATH,
                     sam_index_file  = MM10_UCSC_GENOME,
                     threads         = THREADS,
                     ending_pattern  = 'fq.gz'):
    basename = get_basename(read1, ending_pattern)
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
                        threads           = THREADS,
                        ending_pattern    = 'fq.gz' ):
    basename   = get_basename(read1, ending_pattern)
    hisat2_run = (
          hisat2[
              '-p', threads,
              '--dta',
              '--fr', 
              '-x', hisat2_index_file, 
              '-1', read1, '-2', read2,
          ] | samtools[
              'view',
              '-bSh',
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

def build_BAM_index(input_file, output_file ):
    
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


"""
@aim  this function need sorted bam file. Hence, the BAW or HISAT2 will
      first be carried out and aligned files is save in the working 
      directory.
@params:
      
@return
"""
def run_PICARD_QC_modules( sample_name,
                           ref_genome      = MM10_UCSC_GENOME,
                           ref_flat        = REFFLAT_MM10_UCSC_PICARD,
                           ribo_annotation = RIBO_INTERVAL_LIST_MM10_PICARD,
                           strandness      = STRANDNESS):
    ribo_interval_file = _get_RIBO_file( sample_name, 
                                         ribo_annotation = ribo_annotation)
    _run_picard_CollectRnaSeqMetrics( sample_name,
                                      ribo_interval = ribo_interval_file,
                                      ref_flat      = ref_flat,
                                      strandness    = strandness
                                    )
    _run_picard_CollectAlignmentSummaryMetrics( sample_name,
                                                ref_genome    = MM10_UCSC_GENOME)

    _run_picard_CollectInsertSizeMetrics(sample_name)
    return 1

"""
@aim:  the def is the wrapper for the picard QC module
       CollectRnaSeqMetrics. And perform the QC checking.
@parameters
    sample_name:      get the sample name ( = basename) sample name
                      was setted to be the basename free of ending pattern.
    ribo_inerval:     the interval file which is save for the ribsome location
                      and dynamically genratead by the _ribo_ function.
                      this file is specifically needed by this def.
    ref_flat:         the ref_flat file generated according to the suggestion
                      by BioStar post.
    strandness:       the input parameter transfered from the python script

@return
     none.

"""
def _run_picard_CollectRnaSeqMetrics( sample_name,
                                      ribo_interval,
                                      ref_flat      = REFFLAT_MM10_UCSC_PICARD,
                                      strandness    = STRANDNESS):
    
    picard_run = ( picard[ 'CollectRnaSeqMetrics',
                           'REF_FLAT=', ref_flat,
                           'RIBOSOMAL_INTERVALS=',ribo_interval,
                           'STRAND_SPECIFICITY=',strandness,
                           'CHART_OUTPUT=','null',
                           'METRIC_ACCUMULATION_LEVEL=', 'ALL_READS',
                           'INPUT=',  sample_name + '.bam',
                           'OUTPUT=', sample_name + '.CollectRnaSeqMetrics.picard',
                           'ASSUME_SORTED=','true'] )
    picard_run()
    remove_file = (rm[ribo_interval])
    remove_file()


"""
@aim:  the def is the wrapper for the picard QC module
       CollectAlignmentSummaryMetrics. And perform the 
       corresponding QC checking.
@parameters
    sample_name:      get the sample name ( = basename) sample name
                      was setted to be the basename free of ending pattern.
    ref_genome:       the input parameter which specify the regerence genome
                      abolsote path.

@return
     none.
"""
def _run_picard_CollectAlignmentSummaryMetrics( sample_name,
                                                ref_genome    = MM10_UCSC_GENOME):
    
    picard_run = ( picard[ 'CollectAlignmentSummaryMetrics',
                           'REFERENCE_SEQUENCE=', ref_genome,
                           'INPUT=',  sample_name + '.bam',
                           'OUTPUT=', sample_name + '.CollectAlignmentSummaryMetrics.picard',
                           'EXPECTED_PAIR_ORIENTATIONS=','null']
                  )
    picard_run()

"""
@aim:  the def is the wrapper for the picard QC module
       CollectInsertSizeMetrics. And perform the 
       corresponding QC checking.
@parameters
    sample_name:      get the sample name ( = basename) sample name

@return
     none.
"""

def _run_picard_CollectInsertSizeMetrics(sample_name):
    picard_run = ( picard[ 'CollectInsertSizeMetrics',
                           'INPUT=',  sample_name + '.bam',
                           'OUTPUT=', sample_name + '.CollectInsertSizeMetrics.picard',
                           'HISTOGRAM_FILE=', sample_name + '.CollectInsertSizeMetrics.pdf',
                           'MINIMUM_PCT=', 0.05] )
    picard_run()


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

'''
@aim:
    to  run the fastqc program and return all the checked files.
    this is a inner function and should not be used in the outside
    calling for program usage.
@parameters
    file (String): the raw illumina reads file, single file

@return
    the concanetented string of all files. and create a local
    dir which contains all QC checked results

'''

def _run_FASTQC( file, threads = 1,
                 output_dir = 'fastqc.results'):
    os.makedirs(output_dir, exist_ok = True)
    run_fastqc =  fastqc[ '-o', output_dir,
                           '-f', 'fastq',
                           '-t', threads,
                           file]
    (run_fastqc)() 
    return run_fastqc

'''
@aim:
    to  run the fastqc program and return all the checked files.
    this is a inherited fucntion from the above _run_FASTQC().
@parameters
    file (String/list): the raw illumina reads file, single file or
                        mutilple files in list.

@return
    create a local dir which contains all QC checked results

'''

def run_FASTQC( files, threads = THREADS,
                output_dir = 'fastqc.results'):
    if type(files) is str:
        _run_FASTQC( file, 
                     output_dir = 'fastqc.results')
    elif type(files) is list:
        for file in files:
            _run_FASTQC( file,
                     output_dir = 'fastqc.results')
    else:
        return 'file type error'

'''
@aim:
    to  run the samtools program and generate header file for
    the later ussage. The header file will be combined with 
    the ribo file annnnotion which is used in the picard
    program to determine the Percentage of ribosome file.

@parameters
    file (String): the aligned file in bam foramt, single file
                   and spiece specfic ribo_annotation

@return
    the temp file for picard usage for QC checking.

'''
def _get_RIBO_file( base_name, ribo_annotation = RIBO_INTERVAL_LIST_MM10_PICARD):
    TEMP_FILE      = tempfile.NamedTemporaryFile(dir = os.getcwd())
    TEMP_FILE_NAME = TEMP_FILE.name
    TEMP_FILE.close()
    get_samtool_header = ( samtools[ 'view',
                                     '-H', base_name + '.bam',
                                     '-o', TEMP_FILE_NAME] )
    get_samtool_header()
    
    get_ribo_file  = ( cut[ '-s',
                            '-f', '1,4,5,7,9',
                            ribo_annotation] >> TEMP_FILE_NAME )
    get_ribo_file()
    return TEMP_FILE_NAME

"""
not implemented
"""

def choose_ANNOTATION():
    return 1

