#---
# @author Yisong Zhen
# @since  2018-02-28
# @update 2018-03-19
#---


#---
#
# project aim:
# ChIPseq pipeline in python2 version
# for the sake of MACS2 have to coorespongding python3
# alert!!!
# this is a python2 version, using python2 as the 
# intepreter engine.
# I have to migrate the QC pipeline to python2
# this package or script is designed for usage in 
# I have to use the macs2 env setting in miniconda
# 
#---


#---
#  the script setting path: 
#  nohup python /wa/zhenyisong/sourcecode/core.facility/epigeneticX/cardiac_mef2.ChIP.py &
#  
#---

#---
# Tool: Script to compute the effective genome size: epic-effective
# https://www.biostars.org/p/185948/
#---



"""

raw data is deposited in the 
ArrayExpress
E-MTAB-6213
SRP121284 
@source  Neonatal rat ventricular cardiomyocytes
I wrote an email to Dr. Lehmann and he responded my request.
1: Lehmann LH, Jebessa ZH, Kreusser MM, Horsch A, He T, Kronlage M, Dewenter M,
Sramek V, Oehl U, Krebs-Haupenthal J, von der Lieth AH, Schmidt A, Sun Q,
Ritterhoff J, Finke D, Völkers M, Jungmann A, Sauer SW, Thiel C, Nickel A,
Kohlhaas M, Schäfer M, Sticht C, Maack C, Gretz N, Wagner M, El-Armouche A, Maier
LS, Londoño JEC, Meder B, Freichel M, Gröne HJ, Most P, Müller OJ, Herzig S,
Furlong EEM, Katus HA, Backs J. A proteolytic fragment of histone deacetylase 4
protects the heart from failure by regulating the hexosamine biosynthetic
pathway. Nat Med. 2018 Jan;24(1):62-72. doi: 10.1038/nm.4452. Epub 2017 Dec 11.
PubMed PMID: 29227474.

"""

import os
import sys
import re
import glob2
import tempfile
import metaseq
import numpy as np
from plumbum import local, FG, BG
from plumbum.cmd import cut, rm
from plumbum.commands.processes import ProcessExecutionError, CommandNotFound
from multiprocessing import Pool as multiThreads

THREADS = 6

BWA_INDEX_RN6_PATH  = (
    '/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/BWAIndex/genome.fa' )
RN6_UCSC_GENOME     = (
    '/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa')

'''
https://gist.github.com/danielecook/b5000d45ec2988da9217
how to get the bed file for the low complexity region in genomes
cd 
wget 'http://hgdownload.soe.ucsc.edu/goldenPath/rn6/database/rmsk.txt.gz' -O LCR_rmsk.txt.gz
gunzip -fc LCR_rmsk.txt.gz | grep 'Low_complexity' | cut -f 6,7,8 > low_complexity_region_rn6_rmsk.bed
rm LCR_rmsk.txt.gz
'''
LOW_COMPLEXITY_RN6 = (
      '/wa/zhenyisong/reference/annotation/ChIPQC/low_complexity_region_rn6_rmsk.bed')


def _run_BWA_aligner_SE( read,
                         threads         = THREADS,
                         bwa_index_file  = BWA_INDEX_RN6_PATH,
                         sam_index_file  = RN6_UCSC_GENOME,
                         ending_pattern  = 'fq.gz'):
    try:
        basename = get_basename(read, ending_pattern)
        run_bwa  = ( 
             bwa[ 
                 'mem',
                 '-M',
                 '-t', threads,
                 bwa_index_file,
                 read
             ] | samtools[
                 'view',
                 '-bSh',
                 '-@', THREADS,
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
    
    except ProcessExecutionError:
        print( '''Please check the procedure for sequence
                  alignment run_BWA_aligner module was failed.
               ''')
    except CommandNotFound:
        print('this commnand BWA is not congifured well')
    except:
        print('well, whatever, we failed')
        
    return read

def get_basename(fullname, ending_pattern = 'fq.gz' ):
    basename = os.path.basename(fullname)
    basename = re.sub(ending_pattern,'', basename)
    basename = re.sub('\\.$', '', basename)
    return basename

#---
# Use a Glob() to find files recursively in Python?
#
# https://stackoverflow.com/questions/2186525/
# use-a-glob-to-find-files-recursively-in-python
#
def get_raw_data_names(raw_data_path, ending_pattern = '.fq.gz'):
    
    raw_data_path = raw_data_path + '/**/*' + ending_pattern
    files         = []
    for file in glob2.glob(raw_data_path):
        files.append(file)
    if len(files) == 0:
        raise ValueError('no raw sequencing files found')
    return files


def _run_macs2_peakCalling( treatment, control,
                            bed_filename,
                            file_format = 'BAM',
                            gsize  = 'mm',
                            qvalue = 0.01
                          ):
    try:
       run_macs2  = ( 
                macs2[ 
                      'callpeak',
                      '--treatment', treatment,
                      '--control', control,
                      '--format', file_format,
                      '--gsize', gsize,
                      '--name', bed_filename,
                      '--bdg',
                      '--qvalue', qvalue
                     ] 
                   )
       run_macs2()
    except ProcessExecutionError:
        print( '''Please check the procedure for sequence
                  alignment run_BWA_aligner module was failed.
               ''')
    except CommandNotFound:
        print('this commnand BWA is not configured well')
    except:
        print('well, whatever, we failed') 

def get_BWA_bam_results( chip_seq_files, suffix):
    for file in chip_seq_files:
        _run_BWA_aligner_SE(file, ending_pattern = suffix)

def get_MACS2_results( treats, controls, bed_filename = None,
                       file_format = 'BAM',
                            gsize  = 'mm',
                            qvalue = 0.01 ):
    if bed_filename is None:
        assert len(treats) == len(controls)
        assert len(treats) != 0
        for i in range(len(treats)):
            _run_macs2_peakCalling( treats[i], controls[i],
                                    bed_filename = 
                                       get_base_name(treats[i], '.bam'),
                                    file_format = file_format,
                                    gsize  = gsize,
                                    qvalue = qvalue)
    return None    

    


def split_ChIP_samples(raw_sample_files):
    assert len(raw_sample_files) % 2 == 0
    treat_samples   = raw_sample_files[1::2]
    control_samples = raw_sample_files[::2]
    return treat_samples, control_samples

def get_base_name(fullname, suffix):
    basename = os.path.basename(fullname)
    basename = re.sub(suffix,'', basename)
    basename = re.sub('\\.$','',basename)
    return basename

def set_working_path(working_path):
    os.chdir(working_path)
    return None

"""
how to determine the effective genome size
https://github.com/taoliu/MACS/issues/116
RN6_UCSC_GENOME=
'/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa'
epic-effective --read-length=150 --nb-cpu=4 ${RN6_UCSC_GENOME}
"""
def get_effective_genome_size( genome_file, 
                               threads     = 1, 
                               read_length = 25):
    
    try:
        run_epic_effective = (
              epic_genome_size[
                                '--read-length=', read_length,
                                '--nb-cpu=',threads,
                                genome_file
                              ]                  
                             )
        run_epic_effective()
    except ProcessExecutionError:
        print( '''Please check the procedure for calling epic,
                  get_effective_genome_size module was failed.
               ''')
    except CommandNotFound:
        print('this commnand epic-effective is not configured well')
    except:
        print('well, whatever, we failed') 
    
    return None

def read_MACS2_output():
    return None

def _trim_bam_file( bam_file_name,
                    trim_length):
    sample_name = get_base_name(bam_file_name, suffix = '.bam')
    try:
        run_trim_bam  = (
              trim_bam[
                        'trimBam', 
                        bam_file_name,
                        '-R', trim_length,
                        '--clip'
                      ] | picard[
                        'SortSam', 
                        'INPUT=','/dev/stdin',
                        'OUTPUT=', sample_name + '.trimmed.bam',
                        'SORT_ORDER=','coordinate'
                      ]
             )
        run_trim_bam()
    except ProcessExecutionError:
        print( '''Please check the procedure for calling epic,
                  get_effective_genome_size module was failed.
               ''')
    except CommandNotFound:
        print('this commnand epic-effective is not configured well')
    except:
        print('well, whatever, we failed')
    return sample_name 

def _index_bam_file(filename):
    sample_name = get_base_name(bam_file_name, suffix = '.bam')
    try:
        run_index_bam  = (
              picard[
                      'BuildBamIndex', 
                      'INPUT=', filename,
                      'OUTPUT=', filename + '.bai'
                      ] 
             )
        run_index_bam()
    except ProcessExecutionError:
        print( '''Please check the procedure for calling epic,
                  get_effective_genome_size module was failed.
               ''')
    except CommandNotFound:
        print('this commnand epic-effective is not configured well')
    except:
        print('well, whatever, we failed')
    return sample_name 
'''

I think the following protocol is specific
for python 2.7.
If migrate to python3, the sytax will be changed.
@parameters
    1. bam_files(List)      : the BAM file name list
    2. trim_length(Integer) : the integer to trim the 
                              right bases in the alignment
    3. threads(Integer)     : the threads needed in the
                              setting.
@return
    the sample names(List)
@reference
    python programing for biologists
@update  2018-03-20

'''

def trim_bams_multi_threads( bam_files, 
                             trim_length = 70,
                             threads     = THREADS):
    pool     = multiThreads(threads)
    jobs     = []
    for bam_file in bam_files:
        trim_args = (bam_file, trim_length) 
        job       = pool.apply_async( _trim_bam_file, trim_args ) 
        jobs.append(job)
    results  = []
    for job in jobs:
        result = job.get()
        results.append(result)
    pool.close() 
    pool.join()
    pool.terminate()
    return results


samtools         = local['samtools']
gatk4            = local['gatk-launch']
picard           = local['picard']
fastqc           = local['fastqc']
bwa              = local['bwa']
macs2            = local['macs2']
epic_genome_size = local['epic-effective']
trim_bam         = local['bam']



data_source = '/home/zhenyisong/data/cardiodata/PRJEB23434'
raw_files   = get_raw_data_names(data_source, ending_pattern = 'fastq.gz')
set_working_path(data_source)
get_BWA_bam_results(raw_files, suffix = 'fastq.gz')
bam_files_from_BWA = get_raw_data_names(os.getcwd(), ending_pattern = 'bam')

treats,controls    = split_ChIP_samples(bam_files_from_BWA)


get_MACS2_results( treats, controls, 
                   bed_filename = None,
                   file_format = 'BAM',
                   gsize  = 2.2e9,
                   qvalue = 0.01 )

# metaseq._genomic_signal.supported_formats()


temp = metaseq._genomic_signal.genomic_signal('ERR2194234_summits.bed' ,'bed')
print('now, completed the BWA alignment')