#---
# @author Yisong Zhen
# @since  2018-04-08
# @update 2018-04-13
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
import HTSeq
import collections
from rpy2.robjects.packages import importr 

THREADS = 3
raw_data_path    = '/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs'
analysis_results = '/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs/bwa'
BWA_INDEX_PATH        = (
        '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa' )
REFERENCE_GENOME      = BWA_INDEX_PATH
HG38_UCSC_GTF    = (
        '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf')


def run_BWA_aligner( read1, 
                     read2           = None,
                     library_model   = 'PE',
                     bwa_index_file  = BWA_INDEX_PATH,
                     threads         = THREADS,
                     sorting_method  = 'queryname',
                     output_filename = None,
                     ending_pattern  = 'fq.gz'):
    assert isinstance(read1, str), 'read1 is not string'
    assert isinstance(read2, str) or read2 is None, 'read2 is incorrect value'
    try:
        basename         = get_basename(read1, ending_pattern)
        if output_filename is None:
            output_filename = basename + '.bam'
        if library_model == 'PE' and read2 is not None:
            run_bwa_PE  = ( 
                 bwa[ 
                     'mem',
                     '-M',
                     '-t',threads,
                     bwa_index_file,
                     read1, read2
                 ]  | picard[
                     'SortSam', 
                     'INPUT=','/dev/stdin',
                     'OUTPUT=', output_filename,
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
                 ] | picard[
                     'SortSam', 
                     'INPUT=','/dev/stdin',
                     'OUTPUT=', output_filename,
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
'''

@return a single sample gene read counts in dict type
        the object of collection.Counter class
@update 2018-04-18
'''

def read_gene_counts( bam_file, gene_features_object,
                      library_model  = 'PE',
                      stranded       = 'no',):
    bam_file_object = HTSeq.BAM_Reader( bam_file )
    assert isinstance( gene_features_object, 
                       HTSeq.GenomicArrayOfSets ), 'the object is not right class'
    counts          = collections.Counter( )
    for bundle in HTSeq.pair_SAM_alignments( bam_file_object, bundle = True ):
        if len( bundle ) != 1:
            # Skip multiple alignments
            continue 
        # extract pair
        first_almnt, second_almnt = bundle[0]  
        if ( not hasattr(first_almnt, 'aligned') or 
             not hasattr(second_almnt, 'aligned') ):
            counts[  '_unmapped' ] += 1
            continue
        if not first_almnt.aligned or not second_almnt.aligned:
            counts[  '_unmapped' ] += 1
            continue
        gene_ids = set()
        for iv, val in gene_features_object[ first_almnt.iv ].steps():
            gene_ids |= val
        for iv, val in gene_features_object[ second_almnt.iv ].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[ gene_id ] += 1
        elif len(gene_ids) == 0:
            counts[ '_no_feature' ] += 1
        else:
            counts[ '_ambiguous' ] += 1
    return counts
'''
@parameters
    1. GTF_filename(String)  : the path to GTF file
    2. stranded_info(Binary) : whether the future sam or bam file is
                               stranded or not


'''

def _get_features_file( GTF_filename, stranded_info = False):
    exon_annotation_features   = HTSeq.GenomicArrayOfSets( 'auto', 
                                                stranded = stranded_info )
    gtf_annotation             = HTSeq.GFF_Reader( GTF_filename )
    
    for feature in gtf_annotation:
        if feature.type == 'exon':
            exon_annotation_features[ feature.iv ] += feature.attr['gene_id']
    return exon_annotation_features


def _change_gene_dict_to_table(samples_counts_dict):
    assert isinstance(samples_counts_dict, dict), 'the input is not dictionary'
    samples_list = samples_counts_dict.keys()
    for gene in samples_counts_dict[samples_list[0]].keys():
        



def perform_mRNA_diff_procedure():


temp_file = '/home/zhenyisong/data/results/chenlab/songli/mRNAhumanYs/temp.bam'
#gtf_file  = HTSeq.GFF_Reader( HG38_UCSC_GTF )

features = _get_features_file(HG38_UCSC_GTF)
temp_bam = HTSeq.BAM_Reader( temp_file )
temp = read_gene_counts(temp_file, features)