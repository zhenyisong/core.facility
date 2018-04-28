#---
# @author Yisong Zhen
# @since  2018-04-08
# @update 2018-04-28
#---


# qsub 
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
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
from rpy2.robjects import Formula
import pickle
import cloudpickle

THREADS          = 3
raw_data_path    = '/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs'
analysis_results = '/wa/zhenyisong/results/chenlab/songli/mRNAhumanYs/bwa'
BWA_INDEX_PATH   = (
        '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa' )
REFERENCE_GENOME = BWA_INDEX_PATH
HG38_UCSC_GTF    = (
        '/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf')


def run_BWA_with_limit_memory( read1, 
                               read2           = None,
                               library_model   = 'PE',
                               bwa_index_file  = BWA_INDEX_PATH,
                               threads         = THREADS,
                               sorting_method  = 'coordinate',
                               output_filename = None,
                               ending_pattern  = 'fq.gz'):
    assert isinstance(read1, str), 'read1 is not string'
    assert isinstance(read2, str) or read2 is None, 'read2 is incorrect value'
    sample_name         = get_basename(read1, ending_pattern)
    if output_filename is None:
            output_filename = sample_name + '.bam'
    temp_file_name = _BWA_mapping( read1, 
                                   read2           = read2,
                                   library_model   = library_model,
                                   bwa_index_file  = bwa_index_file,
                                   threads         = threads,
                                   temp_filename   = None)
    _Picard_sorting( temp_file_name, output_filename, sorting_method )
    remove_file = (rm[temp_file_name])
    remove_file()
    return sample_name

def _BWA_mapping( read1, 
                  read2           = None,
                  library_model   = 'PE',
                  bwa_index_file  = BWA_INDEX_PATH,
                  threads         = THREADS,
                  temp_filename   = None):
    try:
        if temp_filename is None:
            TEMP_FILE      = tempfile.NamedTemporaryFile(dir = os.getcwd())
            TEMP_FILE_NAME = TEMP_FILE.name
            temp_filename  = TEMP_FILE_NAME
            TEMP_FILE.close()
        if library_model == 'PE' and read2 is not None:
                run_bwa_PE  = ( 
                     bwa[ 
                         'mem',
                         '-M',
                         '-t',threads,
                         bwa_index_file,
                         read1, read2 
                     ] > temp_filename
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
                 ] > temp_filename
              )
            run_bwa_SE()
    except ProcessExecutionError:
        print( '''Please check the procedure for sequence
                  alignment _BWA_mapping module was failed.
               ''')
        sys.exit(1)
    except CommandNotFound:
        print('this commnand BWA is not congifured well')
        sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        sys.exit(1)
    return temp_filename

def _Picard_sorting( input_file, output_file, sorting_method ):
    try:
        assert ( sorting_method == 'coordinate' or
                 sorting_method == 'queryname'), 'sorting method is wrong'
        sort_mapping_result = ( 
            picard[
                   'SortSam', 
                   'INPUT=', input_file,
                   'OUTPUT=', output_file,
                   'SORT_ORDER=', sorting_method
                  ])
        sort_mapping_result()
    except ProcessExecutionError:
        print( '''Please check the procedure for 
                  _Picard_sorting module was failed.
               ''')
        sys.exit(1)
    except CommandNotFound:
        print('this commnand PICARD is not configured well')
        sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        sys.exit(1)
    return output_file

'''

@return a single sample gene read counts in dict type
        the object of collection.Counter class
@update 2018-04-18
'''

def read_gene_counts( bam_file, gene_features_object,
                      library_model  = 'PE',
                      stranded       = 'no'):
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
    genenames    = []
    gene_sample_count_dic = dict()
    
    for gene in samples_counts_dict[samples_list[0]].keys():
        genenames.append(gene)
        for sample in samples_list:
            gene_sample_count_dic[sample].append(samples_counts_dict[sample]) 
    
    pandas_df = pd.DataFrame( gene_sample_count_dic,
                              index = genenames )
    return pandas_df


def perform_mRNA_diff_procedure( count_dataframe, 
                                 groups,
                                 design):
     sample_groups             = dict()
     categories                = groups
     design                    = Formula('~ design')
     sample_groups['design']   = robjects.IntVector(categories)
     sample_df                 = robjects.DataFrame(sample_groups)

     diff_result = deseq.DESeqDataSetFromMatrix( countData = count_dataframe, 
                                         colData   = sample_df,
                                         design    = design)
     return None
def get_raw_data_names( raw_data_path, 
                        ending_pattern = '.fq.gz'):
    raw_data_path = raw_data_path + '/**/*' + ending_pattern
    files         = []
    try:
        for file in glob.glob(raw_data_path, recursive = True):
            files.append(file)
        if len(files) == 0:
            raise Exception('no raw sequencing files found')
    except Exception as error:
        print(repr(error))
        sys.exit(1)
    finally:
        print('step to get all raw data filenames plus their path')
    return files

def set_working_path(working_dir):
    try:
        if os.path.exists(working_dir) and \
           os.path.isdir(working_dir)  and \
           os.access(working_dir, os.W_OK):
               os.chdir(working_dir)
        elif not os.path.exists(working_dir):
            os.makedirs(working_dir)
            os.chdir(working_dir)
        else:
            raise Exception('cannot create or change the ' + working_dir)
    except Exception as error:
        print(repr(error))
        print('we lose the battle!!!')
        sys.exit(1)
    finally:
        print('change the Python3 script working path' + '\n' + os.getcwd())
    return None
    
def split_PairEnd_files(files):
    assert isinstance(files, list), (
              'the files is not the python list type' )
    files.sort(reverse = True)
    read1_list = files[1::2]
    read2_list = files[::2]
    assert len(read1_list) == len(read2_list), 'the paired-end files is not even number'
    assert len(read2_list) != 0, 'the raw sequencing data have zero files'
    return read1_list, read2_list
      
def get_basename(fullname, ending_pattern = 'fq.gz' ):
    assert  isinstance(fullname, str)
    assert fullname and fullname.strip()
    basename = os.path.basename(fullname)
    basename = re.sub(ending_pattern,'', basename)
    basename = re.sub('\\.$','', basename)
    return basename


picard     = local['picard']
fastqc     = local['fastqc']
bwa        = local['bwa']
samtools   = local['samtools']
    
whole_data_names       = get_raw_data_names( 
                                 raw_data_path, 
                                 ending_pattern = 'clean.fq.gz')
read1_list, read2_list = split_PairEnd_files(whole_data_names)

set_working_path(analysis_results)

'''
for i in range(len(read1_list)):
    run_BWA_with_limit_memory( read1_list[i], 
                               read2_list[i],
                               library_model   = 'PE',
                               bwa_index_file  = BWA_INDEX_PATH,
                               threads         = THREADS,
                               output_filename = None,
                               sorting_method  = 'queryname',
                               ending_pattern  = '.clean.fastq.gz')
'''


features            = _get_features_file(HG38_UCSC_GTF)
sorted_bam_files    = get_raw_data_names( os.getcwd(), 
                                          ending_pattern = '.bam')
whole_sample_counts = dict()
for file in sorted_bam_files:
    sample_name = get_basename( file, ending_pattern = '.bam' )
    whole_sample_counts[sample_name] = read_gene_counts( file,
                                                         features)
#cloudpickle.dumps(whole_sample_counts)
with open('mRNA.songli.pkl', 'wb') as handle:
    pickle.dump( whole_sample_counts, handle, 
                 protocol= pickle.HIGHEST_PROTOCOL)


