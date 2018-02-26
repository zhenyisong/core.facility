#---
#
# project aim:
# ChIPseq pipeline in python3 version
# this package or script is designed for usage in 
# 
#---


#---
#  the script setting path: 
#  python /wa/zhenyisong/sourcecode/polarity.heart/cardiac_motif_ChIP.py
#---

#---
# Tool: Script to compute the effective genome size: epic-effective
# https://www.biostars.org/p/185948/
#---

#---
# @author Yisong Zhen
# @since  2018-02-13
# @update 2018-02-26
#---


#---
# raw data is deposited in the 
# ArrayExpress
# E-MTAB-6213
#---

import os
import sys
import re
import glob
import tempfile
import unittest
from plumbum import local, FG, BG
from plumbum.cmd import cut, rm
from plumbum.commands.processes import ProcessExecutionError, CommandNotFound

THREADS = 6

BWA_INDEX_RN6_PATH  = '/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/BWAIndex/genome.fa'
RN6_UCSC_GENOME     = '/wa/zhenyisong/reference/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa'


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

def get_raw_data_names(raw_data_path, ending_pattern = '.fq.gz'):
    
    raw_data_path = raw_data_path + '/**/*' + ending_pattern
    files         = []
    for file in glob.glob(raw_data_path, recursive = True):
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
    if assertIsNone(bed_filename):
        assert len(treats) == len(controls)
        assert len(treats) != 0
        for i in range(len(treats)):
            _run_macs2_peakCalling( treat1s[i], controls[i],
                                    get_base_name(treats[i]),
                                    file_format = file_format,
                                    gsize  = gsize,
                                    qvalue = qvalue)
    return None    

    


def split_ChIP_samples(raw_sample_files):
    assert len(raw_sample_files) % 2 == 0
    treat_samples   = raw_sample_files[1::2]
    control_samples = raw_sample_files[::2]
    return treat_samples, control_samples

def get_base_name():
    basename = os.path.basename(fullname)
    basename = re.sub(ending_pattern,'', basename)
    basename = re.sub('\\.$','',basename)
    return basename

def set_working_path(working_path):
    os.chdir(working_path)
    return None




samtools = local['samtools']
gatk4    = local['gatk-launch']
picard   = local['picard']
fastqc   = local['fastqc']
bwa      = local['bwa']
macs2    = local['macs2']


data_source = '/home/zhenyisong/data/cardiodata/PRJEB23434'
raw_files   = get_raw_data_names(data_source, ending_pattern = 'fastq.gz')
set_working_path(data_source)
get_BWA_bam_results(raw_files, suffix = 'fastq.gz')
bam_files_from_BWA = get_raw_data_names(os.getcwd(), ending_pattern = 'bam')
treats,controls    = split_ChIP_samples(bam_files_from_BWA)
get_MACS2_results( treats, controls, 
                   bed_filename = None,
                   file_format = 'BAM',
                   gsize  = 'mm',
                   qvalue = 0.01 )

print('now, completed the BWA alignment')