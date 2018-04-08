#---
# @author Yisong Zhen
# @since  2018-04-08
# @update 2018-04-08
#---

#---
# project aim:
# in response the request by sonlgi
# I will write the script to analyze the 
# mRNA differential expressin of her
# human blood mRNA samples
#---

import os
import os.path
import sys
import re
import timeit
from time import sleep
import hashlib
import glob
import psutil
import tempfile
import argparse
import numpy as np
from multiprocessing import Pool as multiThreads
from plumbum import local, FG, BG
from plumbum.cmd import cut, rm, sort, uniq
from plumbum.commands.processes import ProcessExecutionError, CommandNotFound

raw_data_path    = '/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs'
analysis_results = '/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs/bwa'
BWA_INDEX_PATH        = (
        '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa' )
REFERENCE_GENOME      = BWA_INDEX_PATH


def run_BWA_aligner( read1, 
                     read2           = None,
                     library_model   = 'PE',
                     bwa_index_file  = BWA_INDEX_PATH,
                     sam_index_file  = REFERENCE_GENOME,
                     threads         = THREADS,
                     sorting_method  = 'coordinate',
                     middle_name     = None,
                     ending_pattern  = 'fq.gz'):
    assert isinstance(read1, str), 'read1 is not string'
    assert isinstance(read2, str) or read2 is None, 'read2 is incorrect value'
    try:
        basename         = get_basename(read1, ending_pattern)
        output_file_name = None
        if middle_name is None:
            output_file_name = basename + '.bam'
        else:
            output_file_name = basename + '.' + middle_name + '.bam'
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
                     'OUTPUT=', output_file_name,
                     'SORT_ORDER=', sorting_method
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
                     'OUTPUT=', output_file_name,
                     'SORT_ORDER=', sorting_method
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
        #sys.exit(1)
    except CommandNotFound:
        print('this commnand BWA is not congifured well')
        #sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        #sys.exit(1)
    finally:
        print('we have completed BWA alignment module')
        
    return basename
