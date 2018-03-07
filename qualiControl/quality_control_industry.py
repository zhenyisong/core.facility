#---
# project aim:
# quality control pipeline in python3 version
# this package or script is designed for usage in 
# the core facility at GuoZhong.
#---


# the script setting path: 
#  python /wa/zhenyisong/sourcecode/core.facility/qualiControl/quality_control_industry.py
#  
#---



#---
# @author Yisong Zhen
# @since  2018-01-24
# @update 2018-03-07
#---

import os
import sys
import re
import glob
import tempfile
import argparse
import numpy as np
from multiprocessing import Pool as multiThreads
from plumbum import local, FG, BG
from plumbum.cmd import cut, rm
from plumbum.commands.processes import ProcessExecutionError, CommandNotFound


#---
# python style PEP8
# style guide for Python Code:
# constant definition use the UPPER CASE
#---

# claim default setting
#---

STRANDNESS       = None
BWA_INDEX_PATH   = None
REFERENCE_GENOME = None
THREADS          = 3

#--- default CONSTANT end


#---
# annotation and genome files location
# location in linux platform is fixed
# the reference genome sets for various spieces
# were downloaded from iGenomes
#---
MM10_UCSC_GENOME = (
    '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa' )
MM9_UCSC_GENOME  = (
    '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa' )
HG38_UCSC_GENOME = (
    '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa' )
HG19_UCSC_GENOME = (
    '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa' )
RN6_UCSC_GENOME  = (
    '/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa')
MM10_UCSC_GTF    = (
    '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf' )
MM9_UCSC_GTF     = (
    '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf' )
HG38_UCSC_GTF    = (
    '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf')
HG19_UCSC_GTF    = (
    '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf')

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
# 
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
#
# 
# hg38, see more in details mm10
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
# cd /home/zhenyisong/data/reference/annotation/picard
# wget 
# gunzip refFlat.txt.gz; mv refFlat.txt refFlat_hg19.txt
# how to generate hg38_ribosome_interval_list.txt"
#--- 

#---
# now, inlcude mm10, mm9,
#              hg38, hg19,
#              rn6
#---

REFFLAT_MM10_UCSC_PICARD       = (
    '/wa/zhenyisong/reference/annotation/picard/refFlat_mm10.txt' )
REFFLAT_MM9_UCSC_PICARD        = (
    '/wa/zhenyisong/reference/annotation/picard/refFlat_mm9.txt' )
REFFLAT_HG38_UCSC_PICARD       = (
    '/wa/zhenyisong/reference/annotation/picard/refFlat_hg38.txt' )
REFFLAT_HG19_UCSC_PICARD       = (
    '/wa/zhenyisong/reference/annotation/picard/refFlat_hg19.txt' )
REFFLAT_RN6_UCSC_PICARD        = (
    '/wa/zhenyisong/reference/annotation/picard/refFlat_rn6.txt' )
RIBO_INTERVAL_LIST_MM10_PICARD = (
    '/wa/zhenyisong/reference/annotation/picard/mm10_ribosome_interval_list.txt' )
RIBO_INTERVAL_LIST_MM9_PICARD  = (
    '/wa/zhenyisong/reference/annotation/picard/mm9_ribosome_interval_list.txt' )
RIBO_INTERVAL_LIST_HG38_PICARD = (
    '/wa/zhenyisong/reference/annotation/picard/hg38_ribosome_interval_list.txt' )
RIBO_INTERVAL_LIST_HG19_PICARD = (
    '/wa/zhenyisong/reference/annotation/picard/hg19_ribosome_interval_list.txt' )
RIBO_INTERVAL_LIST_RN6_PICARD  = (
    '/wa/zhenyisong/reference/annotation/picard/rn6_ribosome_interval_list.txt' )


#---
# all the following index files were generated 
# by using the script aligner.sh
#
# hisat2 index file location.
# 
# instead,
# BWA index file location, the files were
# bundled with igneome versions.
# and other required annotation files
#---

# HISAT2
HISAT2_INDEX_MM10_PATH = '/wa/zhenyisong/reference/index/mm10'
HISAT2_INDEX_MM9_PATH  = '/wa/zhenyisong/reference/index/mm9'
HISAT2_INDEX_HG38_PATH = '/wa/zhenyisong/reference/index/hg38'
HISAT2_INDEX_HG19_PATH = '/wa/zhenyisong/reference/index/hg19'
HISAT2_INDEX_RN6_PATH  = '/wa/zhenyisong/reference/index/rn6'

# BWA
BWA_INDEX_MM10_PATH    = (
    '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa' )
BWA_INDEX_MM9_PATH     = (
    '/wa/zhenyisong/reference/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa' )
BWA_INDEX_HG38_PATH    = (
    '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa' )
BWA_INDEX_HG19_PATH    = (
    '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa' )
BWA_INDEX_RN6_PATH     = (
    '/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/BWAIndex/genome.fa')

#--- index end




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

#--- plumbum command end

#source deactivate macs2

#multiqc ./

'''
@aim 
    this function use the BWA aligner to genenrate bam files.
    whether the SE or PE model. this BWA aliger wrapper will
    determine the alignment model using the specified parameters.
@parameters needed
   1. threads: (integer) for paraelle computation. This param use 
               the default to THREADS previous defined.
   2. reads1 & read2: (String) the absolute read path in string format
                      file names; two model PE or SE model
                      I have made the read2 is None, for the SE model 
                      choice.  If read1 and read2 both have the input
                      read2 is not None, and library_model choice is
                      PE, whill change the BWA alignment tool to PE
                      alignment procedure.
   3. bwa_index_file (String) : the location of the indexed genome files.
   4. sam_index_file (String) : the required reference genome file used
                                by samtools index procedure. The reference
                                genome file, not the indexed files generated 
                                by samtools or otherway.
   5. ending_pattern (String) : the raw data ending pattern, use of which
                                to extract sample (or base name)
@function
   to perform the BWA alignment and output the
   coordination sorted BAM format alignment files

@test
   run_BWA_aligner(read_1_files[0], read_2_files[0], ending_pattern = '.downsample.fq.gz') 
@return: 
   !!!the read1 and read2 file names (absolute paths)
   !!!in tuple.
   the samlpe_name(or base_name); please refer to get_basename()
@update 03-07-2018

'''

def run_BWA_aligner( read1, 
                     read2 = None,
                     library_model   = 'PE',
                     bwa_index_file  = BWA_INDEX_PATH,
                     sam_index_file  = REFERENCE_GENOME,
                     threads         = THREADS,
                     ending_pattern  = 'fq.gz'):
    try:
        basename = get_basename(read1, ending_pattern)
        if library_model == 'PE' and read2 is not None:
            run_bwa_PE  = ( 
                 bwa[ 
                     'mem',
                     '-M',
                     '-t',threads,
                     bwa_index_file,
                     read1, read2
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
            run_bwa_PE()
        elif library_model == 'SE' and read2 is None:
            run_bwa_SE  = ( 
                 bwa[ 
                     'mem',
                     '-M',
                     '-t',threads,
                     bwa_index_file,
                     read1
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
            run_bwa_SE()
        else:
            raise Exception( """ the reads data input error! 
                                 Maybe the are PE or SE model, 
                                 please confirm the data model and choose the correct
                                 one in the data input!  """)
    
    except ProcessExecutionError:
        print( '''Please check the procedure for sequence
                  alignment run_BWA_aligner module was failed.
               ''')
    except CommandNotFound:
        print('this commnand BWA is not congifured well')
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
    finally:
        print('we have completed BWA alignment module')
        
    return basename

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
@return tuple
   the read1 and read2 file names (absolute paths)
   in tuple.

'''
def run_HISAT2_aligner( read1, read2  = None,
                        library_model = 'PE',
                        hisat2_index_file = HISAT2_INDEX_MM10_PATH,
                        sam_index_file    = MM10_UCSC_GENOME,
                        threads           = THREADS,
                        ending_pattern    = 'fq.gz' ):
    try:
        basename   = get_basename(read1, ending_pattern)
        if library_model == 'PE' and read2 is not None:
            hisat2_run_PE = (
                  hisat2[
                      '-p', threads,
                      '--dta',
                      '--fr', 
                      '-x', hisat2_index_file, 
                      '-1', read1, '-2', read2,
                  ] | samtools[
                      'view',
                      '-bSh',
                      '-@', threads,
                      '-O', 'BAM',
                      '-T', sam_index_file
                  ] | picard[
                      'SortSam', 
                      'INPUT=','/dev/stdin',
                      'OUTPUT=', basename + '.bam',
                      'SORT_ORDER=','coordinate'
                  ]
                )
            hisat2_run_PE()
        elif library_model == 'SE' and read2 is None:
            hisat2_run_SE = (
                  hisat2[
                      '-p', threads,
                      '--dta',
                      '-x', hisat2_index_file, 
                      '-U', read1
                  ] | samtools[
                      'view',
                      '-bSh',
                      '-@', threads,
                      '-O', 'BAM',
                      '-T', sam_index_file
                  ] | picard[
                      'SortSam', 
                      'INPUT=','/dev/stdin',
                      'OUTPUT=', basename + '.bam',
                      'SORT_ORDER=','coordinate'
                  ]
                )
            hisat2_run_SE()
        else:
            raise Exception( """ the reads data input error! 
                                 Maybe the are PE or SE model, 
                                 please confirm the data model 
                                 and choose the correct
                                 one in the data input!  """)
    except ProcessExecutionError:
        print( '''Please check the procedure for sequence
                  alignment run_HISAT2_aligner module was failed.
               ''')
    except CommandNotFound:
        print('this commnand HISAT2 is not configured well')
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
    finally:
        print('we have completed the HISAT2 alignment module')
        
    return basename



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

def _build_BAM_index( sample_name ):
    
    try:
        output_file = sample_name + '.bam' + '.bai'
        run_build_index = ( gatk4[ 
                              'BuildBamIndex', 
                              '--INPUT',  sample_name + '.bam',
                              '--OUTPUT', output_file] )
        run_build_index()
    except ProcessExecutionError:
        print( '''Please check the procedure for index building
                  _build_BAM_index was failed.
               ''')
    except CommandNotFound:
        print('this commnand GATK4 is not configured well')
    except:
        print('well, whatever, we failed')
    return output_file

'''
@aim 
    get the base name of the sequncing raw data for
    interna usage as the sample names.
@parameters
    fullname:         file full name with file path.
    ending_pattern:   the read data ending pattern which 
                      discriminate the raw data from other 
                      non-related files
@return
    the file base-name which may be used as the sample name
    to trace the QC results with the corresponding raw data

'''

def get_basename(fullname, ending_pattern = 'fq.gz' ):
    assert  type(fullname) == str
    assert fullname and fullname.strip()
    basename = os.path.basename(fullname)
    basename = re.sub(ending_pattern,'', basename)
    return basename


"""
@aim  this function need sorted bam file. Hence, the BAW or HISAT2 will
      first be carried out and aligned files is save in the working 
      directory.
@params:
      
@return (String):
      sample name without suffix.
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
    _build_BAM_index(sample_name)
    _run_picard_CollectInsertSizeMetrics(sample_name)
    
    return sample_name

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
    try:
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
    except ProcessExecutionError:
        print( '''Please check the procedure for QC assessment 
                  CollectRnaSeqMetrics module was failed.
               ''')
    except CommandNotFound:
        print('this commnand PICARD is not configured well')
    except:
        print('whatever, we failed in the _run_picard_CollectRnaSeqMetrics')
    finally:
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
    try:
        picard_run = ( picard[ 'CollectAlignmentSummaryMetrics',
                               'REFERENCE_SEQUENCE=', ref_genome,
                               'INPUT=',  sample_name + '.bam',
                               'OUTPUT=', sample_name + '.CollectAlignmentSummaryMetrics.picard',
                               'EXPECTED_PAIR_ORIENTATIONS=','null']
                      )
        picard_run()
    except ProcessExecutionError:
        print( '''Please check the procedure for QC assessment
                  CollectAlignmentSummaryMetrics module was failed.
               ''')
    except CommandNotFound:
        print('this commnand PICARD is not configured well')
    except:
        print('whatever, we failed in the _run_picard_CollectAlignmentSummaryMetrics')


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
    try:
        picard_run = ( picard[ 'CollectInsertSizeMetrics',
                               'INPUT=',  sample_name + '.bam',
                               'OUTPUT=', sample_name + '.CollectInsertSizeMetrics.picard',
                               'HISTOGRAM_FILE=', sample_name + '.CollectInsertSizeMetrics.pdf',
                               'MINIMUM_PCT=', 0.05] )
        picard_run()
    except ProcessExecutionError:
        print( '''Please check the procedure for QC assessment
                  CollectInsertSizeMetrics module was failed.
               ''')
    except CommandNotFound:
        print('this commnand PICARD is not configured well')
    except:
        print('well, whatever, we failed')



'''
@aims
    the function parse the data path and fetch all
    raw data required to perform QC checking
    raw read data are using the relative or absolute
    path.
@parameters
    raw_data_path(string): the raw illumina reads sequencing results
    end_pattern(string)  :   the read data ending pattern which 
                           discriminate them from other non-related
                           files
@return (list)
    the function return the file with file path lists

'''

def get_raw_data_names(raw_data_path, ending_pattern = '.fq.gz'):
    raw_data_path = raw_data_path + '/**/*' + ending_pattern
    files         = []
    try:
        for file in glob.glob(raw_data_path, recursive = True):
            files.append(file)
        if len(files) == 0:
            raise Exception('no raw sequencing files found')
    except Exception as error:
        print(repr(error))
    finally:
        print('step to get all raw data filenames plus their path')
    return files

'''
@aim:
    to get the Pair-end sequencing file list
    read1 and read2 file list.
@parameters
    files: the raw illumina reads file list

@return (list)
    the splitted files in READ1 and READ2 format
    list

'''
def split_PairEnd_files(files):
    read1_list = files[1::2]
    read2_list = files[::2]
    assert len(read1_list) == len(read2_list), 'the paired-end files is not even number'
    assert len(read2_list) != 0, 'the raw sequencing data have zero files'
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
    try:
        run_fastqc =  fastqc[ '-o', output_dir,
                              '-f', 'fastq',
                              '-t', threads,
                              file]
        (run_fastqc)()

    except ProcessExecutionError:
        print( '''Please check the procedure for QC assessment
                  _run_FASTQC was failed.
               ''')
    except CommandNotFound:
        print('this commnand fastqc is not configured well')
    except:
        print('well, whatever, we failed')
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

@depracate? I will use the multi-thread model to check the raw
            data QC.

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
    to  run the fastqc program and use the multiple threads.
    this is a inherited fucntion from the above _run_FASTQC()
    and redesigned and will take place of the old way, single
    thread way.
@parameters
    file (String/list): the raw illumina reads file, single file or
                        mutilple files in list.

@return
    create a local dir which contains all QC checked results

'''

def run_multiThreads_FASTQC( files, threads = THREADS,
                             output_dir     = 'fastqc.results'):
    
    if type(files) is str:
        _run_FASTQC( file, 
                     output_dir = 'fastqc.results')
    elif type(files) is list:
        pool        = multiThreads(THREADS)
        threads     = np.repeat(np.array([1]), [len(files)])
        output      = np.repeat(np.array([output_dir]), [len(files)])
        fastqc_args = zip(files, threads,output)
        pool.starmap( _run_FASTQC, fastqc_args ) 
        pool.close() 
        pool.join()
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

@return (String)
    the temp file name for picard usage for QC checking.

'''

def _get_RIBO_file( base_name, ribo_annotation = RIBO_INTERVAL_LIST_MM10_PICARD):
    TEMP_FILE      = tempfile.NamedTemporaryFile(dir = os.getcwd())
    TEMP_FILE_NAME = TEMP_FILE.name
    TEMP_FILE.close()
    
    try:
        get_samtool_header = ( samtools[ 'view',
                                         '-H', base_name + '.bam',
                                         '-o', TEMP_FILE_NAME] )
        get_samtool_header()
        
        get_ribo_file  = ( cut[ '-s',
                                '-f', '1,4,5,7,9',
                                ribo_annotation] >> TEMP_FILE_NAME )
        get_ribo_file()
        
    except ProcessExecutionError:
        print( '''Please check the procedure to create RIBO annotation
                  _get_RIBO_file was failed.
               ''')
    except CommandNotFound:
        print('this commnand samtools or unix shell is not configured well')
    except:
        print('well, whatever, we failed')
    return TEMP_FILE_NAME

"""
not implemented
"""

def choose_ANNOTATION():
    return 1

"""
these are mimic switch phrase
see here the original post about this
issue.
at stackoverflow
"""
def switch_genome_build(choice):
    
    return {
        'hg38': HG38_UCSC_GENOME,
        'hg19': None,
        'mm10': MM10_UCSC_GENOME,
        'rn6' : RN6_UCSC_GENOME
    }[choice]

def switch_BWA_index(choice):
    
    return {
        'hg38': BWA_INDEX_HG38_PATH,
        'hg19': None,
        'mm10': BWA_INDEX_MM10_PATH,
        'rn6' : BWA_INDEX_RN6_PATH
    }[choice]

def switch_strandness(choice):
    
    return {
        'NONE': 'NONE',
        'FR'  : 'FIRST_READ_TRANSCRIPTION_STRAND',
        'RF'  : 'SECOND_READ_TRANSCRIPTION_STRAND'
    }[choice]

def set_working_path(working_dir):
    try:
        if os.path.exists(working_dir) and \
           os.path.isdir(working_dir)  and \
           os.access(working_dir, os.W_OK):
               os.chdir(working_dir)
        else:
            raise Exception('canot create or change the working_dir')
    except Exception as error:
        print(repr(error))
    finally:
        print('change the working path')
    return None

def perform_mRNA_QCtask(params_object):
    param_dict = vars(param_parser.parse_args())
    
    #---
    # now get all required parameters from script inputs
    #---
    
    STRANDNESS       = switch_strandness( param_dict['strandness'] )
    BWA_INDEX_PATH   = switch_BWA_index( param_dict['genome_build'] )
    REFERENCE_GENOME = switch_genome_build( param_dict['genome_build'] )
    file_suffix      = param_dict['name_pattern']
    QC_type          = param_dict['QC_type']
    data_path        = param_dict['data_path']
    library_model    = param_dict['library_model']
    THREADS          = param_dict['threads']
    working_path     = param_dict['working_path']

    whole_data_names = get_raw_data_names( 
                                 os.getcwd(), 
                                 ending_pattern = file_suffix)
    set_working_path(working_path)
    run_FASTQC(whole_data_names)
    if library_model == 'PE':
        read1_list, read2_list = split_PairEnd_files(whole_data_names)
        for i in range(len(read1_list)):
            sample_name = run_BWA_aligner( read1_list[i], 
                                           read2_list[i],
                                           ending_pattern = file_suffix)
            run_PICARD_QC_modules(sample_name)
    elif library_model == 'SE':
        for i in range(len(whole_data_names)):
            sample_name = run_BWA_aligner( whole_data_names, 
                                           ending_pattern = file_suffix)
            run_PICARD_QC_modules(sample_name)
    return None


"""
@aim   open the debugg model to find and test function
       now the script get the input parameters
       which are specified by the end user.
@param None

@return (list) paired end seq files

"""
def debug_model(ending_pattern = '.downsample.fq.gz'):
    
    raw_data_pattern = '/home/zhenyisong/data/temp/test'
    working_dir      = '/home/zhenyisong/data/temp/test'
    os.chdir(working_dir)

    raw_data_PEfile_list = get_raw_data_names( 
                                 os.getcwd(), 
                                 ending_pattern = ending_pattern)
    run_FASTQC(data_PEfile_list)
    split_PairEnd_files(data_PEfile_list)
    return  data_PEfile_list 



#---
# how to get the outside parameter from script 
# command line
# simple argparse example wanted: 1 argument, 3 results
# https://stackoverflow.com/questions/7427101/simple-argparse-example-wanted-1-argument-3-results
# in the following setting
# I created 6 parameters which can be recieved from the srcipt inputs
#---


param_parser = argparse.ArgumentParser()
param_parser.add_argument( '-q', '--QC-type', default = 'mRNA', 
                           choices = [ 'mRNA', 'miRNA', 'lncRNA',
                                       'ChIPseq','DNAseq'])

param_parser.add_argument( '-s', '--strandness', default = 'NONE',
                           choices = [ 'NONE', 'FR','RF'],
                           help = """ this will set the strandness 
                                      in mNRA or lncRNA library
                                      construction method, which 
                                      is valuable to mapping and QC strategy """ )              
param_parser.add_argument( '-d', '--data-path', default = "./", required = False,
                           help = """ this will set the raw data path, 
                                      which user can
                                      read sequencing data from """ )
param_parser.add_argument( '-g', '--genome-build', default = 'mm10', 
                           choices = [ 'hg38', 'hg19', 'mm10',
                                       'mm9','rn6'],
                           help = """ this will set the genome version to quality, 
                                      control the data """ )
param_parser.add_argument( '-l', '--library-model', default = 'PE', 
                           choices = [ 'PE','SE','MA'],
                           help = """ this will set the read whether is paired or not, 
                                      single end - SE, paired end - PE """ )
param_parser.add_argument( '-n', '--name-pattern', default = '.fq.gz', required = True,
                           help = """ this will set the read files suffix """ )
param_parser.add_argument( '-w', '--working-path', default = './', required = False,
                           help = """ this will set the output directory """ )
param_parser.add_argument( '-t', '--threads', default = 3, type = int, required = False,
                           help = """ this will set the thread number which is used in 
                                      fastqc module and bwa module """ )
param_parser.add_argument( '-a', '--aligner', default = 'BWA', choices = ['BWA','HISAT2'],
                           help = """ this will set the alignment aloroithm is used in 
                                      alignment procedure when to genenrate BAM files """ )



print (param_dict)

sys.exit()




#---
# to test this script is good,
# I downsample the data to 1% of the raw data
# source activate biotools
# cp /home/zhenyisong/data/results/chenlab/xiaoning/data/A_1/*.fq.gz /home/zhenyisong/data/temp/test
# seqtk sample -s 100 A_1_R1.fq.gz 0.01 | gzip - > A_1_R1.downsample.fq.gz
# seqtk sample -s 100 A_1_R2.fq.gz 0.01 | gzip - > A_1_R2.downsample.fq.gz
# this  will speed up the developemnt of the pyton3 QC pipeline
#---

data_PEfile_list   = debug_model()
run_FASTQC(data_PEfile_list)
read1_list, read2_list   = split_PairEnd_files(data_PEfile_list )
run_BWA_aligner(data_PEfile_list[0], data_PEfile_list[1], ending_pattern = '.downsample.fq.gz')
run_HISAT2_aligner(data_PEfile_list[0], data_PEfile_list[1], ending_pattern = '.downsample.fq.gz')

run_PICARD_QC_modules('A_1_R1')