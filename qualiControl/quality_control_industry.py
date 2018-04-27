#---
# @author Yisong Zhen
# @since  2018-01-24
# @update 2018-04-27
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
python  /wa/zhenyisong/sourcecode/core.facility/qualiControl/quality_control_industry.py \
-n '.downsample.fq.gz' -g 'mm10' -l 'PE' -s 'NONE' -w /home/zhenyisong/data/temp \
-d /home/zhenyisong/data/cardiodata/test/SRP109298

clean the all the output

rm -rf *.bam *.bai *.picard *.pdf fastqc*
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
# since 2018-03-09, I prefer use the PICARD module to
# take place of the RSeQC due to the reason which the RSeQC is
# maitained using the Python2 library. I have to discard this package.
# I have to use the Python3 to develop and maitain the current QC pipeline
# Now the annotation file updation is slowed anf incomplete.
# @since   2018-03-09
# @update  2018-03-09  
#---

RRNA_HG38_RSEQC                = (
     '/wa/zhenyisong/reference/annotation/RSeQC/hg38_rRNA.bed' )
HOUSE_KEEPING_GENES_HG38_RSEQC = (
    '/wa/zhenyisong/reference/annotation/RSeQC/hg38.HouseKeepingGenes.bed' )
BASIC_GENES_GENCODE_HG38_RSEQC = (
    '/wa/zhenyisong/reference/annotation/RSeQC/hg38_GENCODE_v24_basic.bed' )
RRNA_MM10_RSEQC                = (
    '/wa/zhenyisong/reference/annotation/RSeQC/mm10_rRNA.bed' )
HOUSE_KEEPING_GENES_MM10_RSEQC = (
    '/wa/zhenyisong/reference/annotation/RSeQC/mm10.HouseKeepingGenes.bed' )
BASIC_GENES_GENCODE_MM10_RSEQC = (
    '/wa/zhenyisong/reference/annotation/RSeQC/mm10_GENCODE_VM11_basic.bed' )


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
# @update  2018-03-09
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


'''
@reference
    1. http://genomespot.blogspot.com/2016/06/screen-for-mycoplasma-contamination-in.html
    2. https://www.ncbi.nlm.nih.gov/pubmed/25712092
       1: Olarerin-George AO, Hogenesch JB. Assessing the prevalence of mycoplasma
       contamination in cell culture via a survey of NCBI's RNA-seq archive. Nucleic
       Acids Res. 2015 Mar 11;43(5):2535-42.
    3. mycoplasma_genomes
       download the genomes using the blog links(1, with release 38 updation)

#Acholeplasma laidlawii PG-8A (NC 010163.1).
wget ftp://ftp.ensemblgenomes.org/pub/\
release-38/bacteria/fasta/bacteria_14_collection/\
acholeplasma_laidlawii_pg_8a/dna/\
Acholeplasma_laidlawii_pg_8a.ASM1878v1.dna.toplevel.fa.gz
#Mycoplasma fermentans M64 (NC 014921.1)
#Mycoplasma hominis ATCC23114 (NC 013511.1)
#M. hyorhinisMCLD(NC 017519.1)

echo '>Marginini' > Marginini.fa;
zcat Mycoplasma_arginini_7264.version_1.0.dna.toplevel.fa.gz \
| grep -v '>' >> Marginini.fa

echo '>Mhyorhinis' > Mhyorinis.fa;
zcat Mycoplasma_hyorhinis_sk76.ASM31363v1.dna.toplevel.fa.gz \
| grep -v '>' >> Mhyorinis.fa

echo '>Alaidlawii' > Alaidlawii.fa;
zcat Acholeplasma_laidlawii_pg_8a.ASM1878v1.dna.toplevel.fa.gz \
| grep -v '>' >> Alaidlawii.fa

echo '>Mfermentans' > Mfermentans.fa;
zcat  Mycoplasma_fermentans_pg18.ASM20973v1.dna.toplevel.fa.gz \
| grep -v '>' >> Mfermentans.fa

echo '>Mhominis' > Mhominis.fa;
zcat Mycoplasma_hominis_atcc_23114.ASM8586v1.dna.toplevel.fa.gz \
| grep -v '>' >> Mhominis.fa

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000420105.1_ASM42010v1/GCF_000420105.1_ASM42010v1_genomic.fna.gz
echo '>Morale' > Morale.fa
zcat GCF_000420105.1_ASM42010v1_genomic.fna.gz \
|grep -v '>' >> Morale.fa


echo '>Morale' > Morale.fa;
zcat GCF_000420105.1_ASM42010v1_genomic.fna.gz \
|grep -v '>' >> Morale.fa


rm Myco.fa 2>/dev/null;
cat *fa > Myco.fa;
seqret Myco.fa myco.fa
mv myco.fa Myco.fa
for i in *fa ; do
  bwa index $i
done

'''
MYCOPLASMA_GENOMES           = '/wa/zhenyisong/reference/mycoplasma_genomes/Myco.fa'
MYCOPLASMA_GENOMES_BWA_INDEX = MYCOPLASMA_GENOMES


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
#        the linux shell commands,
#     -- define these linux commands in shell way
#---

hisat2     = local['hisat2']
samtools   = local['samtools']
gatk4      = local['gatk-launch']
picard     = local['picard']
fastqc     = local['fastqc']
bwa        = local['bwa']
bam2fastq  = local['bamToFastq']

#--- plumbum command end

#source deactivate macs2

#multiqc ./

'''
@aim 
    this function use the BWA aligner to genenrate bam files.
    whether the SE or PE model. this BWA aliger wrapper will
    determine the alignment model using the specified parameters.
@parameters needed
  
   1. reads1 & read2: (String) the absolute read path in string format
                      file names; two model PE or SE model
                      I have made the read2 is None, for the SE model 
                      choice.  If read1 and read2 both have the input
                      read2 is not None, and library_model choice is
                      PE, whill change the BWA alignment tool to PE
                      alignment procedure.
   2. library_model (string)  : PE or SE model. If other is specified,
                                the we should modify this function.
   3. bwa_index_file (String) : the location of the indexed genome file by
                                BWA alorithm which should be in line with 
                                with the BWA algorithm version.!!
   4. sam_index_file (String) : the required reference genome file used
                                by samtools index procedure. The reference
                                genome file, not the indexed files generated 
                                by samtools or otherway.
   5. threads: (integer) for paraelle computation. This param use 
               the default to THREADS previous defined.
     
   6. ending_pattern (String) : the raw data ending pattern, use of which
                                to extract sample (or base name)
@function
   to perform the BWA alignment and output the
   coordination sorted BAM format alignment files

@test
   run_BWA_aligner(read_1_files[0], read_2_files[0], ending_pattern = '.downsample.fq.gz') 
@return (String): 
   !!!the read1 and read2 file names (absolute paths)
   !!!in tuple.
   the samlpe_name(or base_name); please refer to get_basename();
   and BWA aligner generated file ( with .bam suffix) will be saved 
   in the current working directory.
@update 04-17-2018

'''

def run_BWA_aligner( read1, 
                     read2           = None,
                     library_model   = 'PE',
                     bwa_index_file  = BWA_INDEX_PATH,
                     threads         = THREADS,
                     sorting_method  = 'coordinate',
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


'''
@parameters needed
   1. read1 (String) : the file name with the path.
   2. read2 (String) : default is None. If PE model is
                       set, also accept the file name 
                       for raw read data with path or 
                       nake file name
   3. library_model(String) : PE or SE or other choice.
                              if other choice, mush re-implement
                              the function procedure.
   4. hisat2_index_file (String) : index file name with path.
   5. sam_index_file (String)    : reference genome file which will
                                   be used to gnerate BAM file
   6. threads (Integer)          : thread number for high performance
   7. ending_pattern             : the raw read data suffix
@function
   to perform the HISAT2 alignment and output the
   coordination sorted BAM format alignment files

@test
   run_HISAT2_aligner(read_1_files[0], read_2_files[0]) 
@return String
   the sample_name, only the sample name without any 
   path information.
   and all BAM files will be dumped into the current working 
   directory.
@update  2018-03-14

'''
def run_HISAT2_aligner( read1, read2  = None,
                        library_model = 'PE',
                        hisat2_index_file = HISAT2_INDEX_MM10_PATH,
                        sam_index_file    = MM10_UCSC_GENOME,
                        threads           = THREADS,
                        ending_pattern    = 'fq.gz' ):
    
    assert isinstance(read1, str)
    assert isinstance(read2, str) or read2 is None
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
        sys.exit(1)
    except CommandNotFound:
        print('this commnand HISAT2 is not configured well')
        sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        sys.exit(1)
    finally:
        print('we have completed the HISAT2 alignment module')
        
    return basename

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
        #sys.exit(1)
    except CommandNotFound:
        print('this commnand BWA is not congifured well')
        #sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        #sys.exit(1)
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
        #sys.exit(1)
    except CommandNotFound:
        print('this commnand PICARD is not configured well')
        #sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        #sys.exit(1)
    return output_file


'''
@aim
    using the mapped file (bam) to extract the fastq seqeuncing
    and mapping against another genome.
@parameters:
    1. filename(String)         : the bam file name
    2. bwa_index_file(string)   : the index file path for the BWA
    3. output_file_name(string) :
    4. threads(Integer)         : the thread number
    5. suffix(string)           : the file suffix so
                                  to extract sample name
@return
    the output_file_name(string)
    and save the bam statistic flie in the output.

@update 2018-03-22
'''
def _run_BWA_reversed_mapping( bam_filename,
                               suffix           = '.bam',
                               bwa_index_file   = BWA_INDEX_PATH,
                               output_file_name = None,
                               threads          = THREADS):
    sample_name = get_basename(bam_filename, ending_pattern = suffix)
    assert isinstance(sample_name, str)
    assert sample_name and sample_name.strip()
    if output_file_name is None:
        output_file_name = sample_name + '.rev.bam'
    try:
        run_re_mapping = ( 
                bam2fastq[
                    '-i', bam_filename,
                    '-fq', '/dev/stdout'
                ] | bwa[
                    'mem',
                    '-t', threads,
                    bwa_index_file,
                    '-'
                ] | samtools[
                    'view',
                    '-F', 4,
                    '-uSh',
                    '-'
                ] | samtools[
                    'sort',
                    '-o', output_file_name,
                    '-'
                ]
        )
        run_re_mapping()
    except ProcessExecutionError:
        print( '''Please check the procedure _run_BWA_reversed_mapping
                  module was failed.
               ''')
        sys.exit(1)
    except CommandNotFound:
        print('this commnand _run_BWA_reversed_mapping is not configured well')
        sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        sys.exit(1)
    finally:
        print('we have completed the _run_BWA_reversed_mapping alignment module')
    return output_file_name
'''
@aim 
   to compute the statistics of the bam file
@parameters
    1.bam_file(string)  :  a single bam file with path
    2.suffix(String)    : some specific ending pattern
    3.output_file_name  : the statistic output flat file
@return(String):
    the saved flat file for the statistics
@update  2018-03-26
'''

def _extract_samtool_stats( bam_file,
                            suffix = 'myco.bam',
                            output_file_name = None):
    try:
        sample_name = get_basename(bam_file, ending_pattern = suffix)
        if output_file_name is None:
            output_file_name = sample_name + '.' + suffix + '.stats.txt'
        extract_bam_stats = (
                    samtools[
                        'view',
                        '-q', 20,
                        bam_file
                    ] | cut[
                        '-f', 3
                    ] | sort | uniq['-c'] > output_file_name )
        extract_bam_stats()
    except ProcessExecutionError:
        print( '''Please check the procedure for sequence
                  _extract_samtool_stats was failed.
               ''')
        sys.exit(1)
    except CommandNotFound:
        print('this commnand unix shells are not configured well')
        sys.exit(1)
    except Exception as error:
        print('Caught this error, we falied: ' + repr(error))
        sys.exit(1)
    finally:
        print('we have completed the _extract_samtool_stats module')
    return output_file_name

'''
@aim 
    get the bam index used the GATK4 program
    this will speed the program scanning of the 
    bam result.
@parameters
    sample_name (Stirng): samle name, which is in line with
                          the corresponding bam file prefix.
                          and bam file should be located in the
                          same working directory.
@return (String)
    bam indexed file name. the returned name will be
    samplename + '.bam' + '.bai'
    the indexed bam file (.bai) is generated and saved in
    the working directory.
@update   2018-03-09

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
    finally:
        print('now the hiddden procedure to create the bai index file is completed')
    return output_file

'''
@aim 
    get the base name of the sequncing raw data for
    interna usage as the sample names. the basename is stripped 
    of the '.' if it is at the end.
@parameters
    fullname(string)        : file full name with file path.
    ending_pattern(String)  : the read data ending pattern which 
                              discriminate the raw data from other 
                              non-related files
@return
    the file base-name which may be used as the sample name
    to trace the QC results with the corresponding raw data
@update   2018-03-14

'''

def get_basename(fullname, ending_pattern = 'fq.gz' ):
    assert  isinstance(fullname, str)
    assert fullname and fullname.strip()
    basename = os.path.basename(fullname)
    basename = re.sub(ending_pattern,'', basename)
    basename = re.sub('\\.$','', basename)
    return basename


"""
@aim  this function need sorted bam file. Hence, the BAW or HISAT2 will
      first be carried out and aligned files is save in the working 
      directory.
@params:
      1. sample_name(String): the naked file namd devoid of path and suffix
      2. ref_genome(String) : the reference genome (filename) used by PICARD
                              CollectAlignmentSummaryMetrics module.
      3. ref_flat(String)   : the file name for refFlat file used by PICARD.
      4. ribo_annotation(String)    : the ribbon annotation filename used by PICARD.
      5. strandness(String)         : the FR or RF or NONE used by PICARD to
                                      infer the library strandness.
@return (String):
      sample name without suffix. 
      and all QC metrics will be save with suffix
      (.picard) in the current working directory.
@update 2018-03-09
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
                                                ref_genome    = ref_genome)
    _build_BAM_index(sample_name)
    _run_picard_CollectInsertSizeMetrics(sample_name)
    
    return sample_name

"""
@aim:  the def is the wrapper for the picard QC module
       CollectRnaSeqMetrics. And perform the QC checking.
@parameters
    1. sample_name(String) : get the sample name ( = basename) sample name
                            was setted to be the basename free of ending pattern.
    2. ribo_inerval(String):  the interval file which is save for the ribsome location
                              and dynamically genratead by the _ribo_ function.
                              this file is specifically needed by this def.
    3. ref_flat(String)    :  the ref_flat file generated according to the suggestion
                              by BioStar post. <file name for ref_flat.
    4. strandness(String)  :  the input parameter transfered from the python script

@return
     sample_name(string).
@update  2018-03-09

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
    return sample_name



"""
@aim:  the def is the wrapper for the picard QC module
       CollectAlignmentSummaryMetrics. And perform the 
       corresponding QC checking.
@parameters
    1. sample_name:   get the sample name ( = basename) sample name
                      was setted to be the basename free of ending pattern.
    2. ref_genome:    the input parameter which specify the regerence genome
                      abolsote path.

@return
     sample_name(string)
@update 2018-03-09
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
        sys.exit(1)
    except CommandNotFound:
        print('this commnand PICARD is not configured well')
        sys.exit(1)
    except:
        print('whatever, we failed in the _run_picard_CollectAlignmentSummaryMetrics')
        sys.exit(1)
    finally:
        print('we have completed CollectAlignmentSummaryMetrics module')
    return sample_name


"""
@aim:  the def is the wrapper for the picard QC module
       CollectInsertSizeMetrics. And perform the 
       corresponding QC checking.
@parameters
    1. sample_name:   input is the sample name ( = basename) sample name
                      and all required bam file should be in the same
                      working directory.

@return (string)
     sample_name
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
        sys.exit(1)
    except CommandNotFound:
        print('this commnand PICARD is not configured well')
        sys.exit(1)
    except:
        print('well, whatever, we failed')
        sys.exit(1)
    finally:
        print('we have completed the CollectInsertSizeMetrics module')
    return sample_name


'''
@aims
    the function parse the data path and fetch all
    raw data required to perform QC checking
    raw read data are using the relative or absolute
    path.
@parameters
    1. raw_data_path(string): the raw illumina reads sequencing results
                              the file path, not the file name.
    2. end_pattern(string)  : the read data ending pattern which 
                              discriminate them from other non-related
                              files
@return (list)
    the function return the file with file path lists
@update  2018-03-09

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
        sys.exit(1)
    finally:
        print('step to get all raw data filenames plus their path')
    return files

'''
@aim:
    to get the Pair-end sequencing file list
    read1 and read2 file list. This function
    seperate the PE reads files.
@parameters
    1. files: the raw illumina paired-end reads file list

@return (list)
    the splitted files in READ1 and READ2 format
    list. Two list contain read1 and read2 seperatedly.
@update 2018-03-14

'''
def split_PairEnd_files(files):
    assert isinstance(files, list), (
              'the files is not the python list type' )
    files.sort(reverse = True)
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
    1. file (String)      : the raw illumina reads file, single file
                            filename with path.
    2. threads (Interger) : I enforced the thread to be 1.
    3. output_dir(String) : the fastqc result output directory.
                            this will be the sub-dir in the current
                            working directory.

@return
    the single read file used by fastqc program. 
    and create a local
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
        sys.exit(1)
    except CommandNotFound:
        print('this commnand fastqc is not configured well')
        sys.exit(1)
    except:
        print('well, whatever, we failed')
        sys.exit(1)
    return file

'''
@aim:
    to  run the fastqc program and return all the checked files.
    this is a inherited fucntion from the above _run_FASTQC().
@parameters
    1. file (String/list): the raw illumina reads file, single file or
                        mutilple files in list.
    2. output_dir        : fastqc output dir.

@return (None)
    create a local dir which contains all QC checked results
@update   2018-03-09

@depracate? I will use the multi-thread model to check the raw
            data QC.


'''

def run_FASTQC( files, threads = THREADS,
                output_dir = 'fastqc.results'):
    try:
        if type(files) is str:
            _run_FASTQC( file, 
                         output_dir = 'fastqc.results')
        elif type(files) is list:
            for file in files:
                _run_FASTQC( file,
                         output_dir = 'fastqc.results')
        else:
            raise Exception('file type error in the fastqc program input')
    except Exception as error:
        print(repr(error))
    finally:
        print('we have completed the FASTQC procedure')
    return None


'''
@aim:
    to  run the fastqc program and use the multiple threads.
    this is a inherited fucntion from the above _run_FASTQC()
    and redesigned and will take place of the old way, single
    thread way.
@parameters
    1. file (String/list    : the raw illumina reads file, single file or
                              mutilple files in list.
    2. threads(Integer)     : the threads needed for fastqc program
    3. output_dir(String)   : the factqc output result position

@return (None)
    None.
    and create a local dir which contains all QC checked results
@update  2018-03-14

'''

def run_multiThreads_FASTQC( files, threads = THREADS,
                             output_dir     = 'fastqc.results'):
    try:
        if type(files) is str:
            _run_FASTQC( file, 
                         output_dir = output_dir)
        elif type(files) is list:
            pool        = multiThreads(THREADS)
            threads     = np.repeat(np.array([1]), [len(files)])
            output      = np.repeat(np.array([output_dir]), [len(files)])
            fastqc_args = zip(files, threads, output)
            pool.starmap( _run_FASTQC, fastqc_args ) 
            pool.close() 
            pool.join()
        else:
            raise Exception('file type error in the run_multiThreads_FASTQC module')
    except Exception as error:
        print('multi-threads-fastqc running error' + repr(error))
        sys.exit(1)
    finally:
        print('we have completed the multi-threads-fastqc module')
    return None



'''
@aim:
    to  run the samtools program and generate header file for
    the later ussage. The header file will be combined with 
    the ribo file annnnotion which is used in the picard
    program to determine the Percentage of ribosome file.

@parameters
    1. base_name(String): or sample_name 
                          the aligned file in bam foramt, 
                          single file
                          and spiece specfic ribo_annotation
    2. ribo_annotation(String) : the predefined constant file name 
                                 with path.

@return (String)
    the temp file name for picard usage for QC checking.
    and a temporary file saved for robiosome annotation is created on-the-fly.
@update 2018-03-09

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
        sys.exit(1)
    except CommandNotFound:
        print('this commnand samtools or unix shell is not configured well')
        sys.exit(1)
    except:
        print('well, whatever, we failed _get_RIBO_file')
        sys.exit(1)
    
    return TEMP_FILE_NAME

'''
@aim to check the contamination of mycoplasma
        in the cell lines.
        see this paper:
        Assessing the prevalence of mycoplasma contamination 
        in cell culture via a survey of NCBI's RNA-seq archive
        PMID: 25712092
        see the experimental log. see the original blog
        and his script.
        http://genomespot.blogspot.jp/2016/06/screen-for-mycoplasma-contamination-in.html

@parameters:
    filenames(String): only accept single end file;
                       if paired end, then the one read file is ok!
                       But have to input splitted file list;
@return
    None
@update  2018-04-25

@replication 
from the Blog script, minor modification!!

cd /home/zhenyisong/data/cardiodata/mycoplasma
MYCOPLASMA_GENOMES='/wa/zhenyisong/reference/mycoplasma_genomes/Myco.fa'
MYCOPLASMA_GENOMES_BWA_INDEX=${MYCOPLASMA_GENOMES}

bwa mem -M -t 4  ${MYCOPLASMA_GENOMES_BWA_INDEX} SRR488569.fastq.gz > temp.sam
samtools view -bSh -@ 4 -O BAM -T ${MYCOPLASMA_GENOMES} -o temp.bam temp.sam
picard SortSam INPUT=temp.bam OUTPUT=temp.sorted.bam SO='coordinate'
picard BuildBamIndex I=temp.sorted.bam 
samtools view -q 20 temp.sorted.bam | cut -f 3 | sort | uniq -c > temp.stats

os.chdir('/home/zhenyisong/data/cardiodata/mycoplasma')
REFERENCE_GENOME = HG38_UCSC_GENOME
BWA_INDEX_PATH   = BWA_INDEX_HG38_PATH
THREADS = 3

check_mycoplasma_contamination( 'SRR944282.fastq.gz',
                                 read2           = None,
                                 library_model   = 'SE',
                                 myco_bwa_index_file  = MYCOPLASMA_GENOMES_BWA_INDEX,
                                 myco_sam_index_file  = MYCOPLASMA_GENOMES,
                                 bwa_index_file       = BWA_INDEX_PATH,
                                 sam_index_file       = REFERENCE_GENOME,
                                 threads         = THREADS,
                                 sorting_method  = 'coordinate',
                                 ending_pattern  = 'fastq.gz' )
                                  

'''

def check_mycoplasma_contamination( read1,
                                    read2           = None,
                                    library_model   = 'SE',
                                    myco_bwa_index_file  = MYCOPLASMA_GENOMES_BWA_INDEX,
                                    myco_sam_index_file  = MYCOPLASMA_GENOMES,
                                    bwa_index_file       = BWA_INDEX_PATH,
                                    sam_index_file       = REFERENCE_GENOME,
                                    threads         = THREADS,
                                    sorting_method  = 'coordinate',
                                    ending_pattern  = 'fq.gz' ):
    sample_name      = get_basename(read1)
    output_file_name = sample_name + '.myco.bam'
    run_BWA_with_limit_memory( read1, 
                               read2           = None,
                               library_model   = 'SE',
                               bwa_index_file  = myco_bwa_index_file,
                               threads         = threads,
                               sorting_method  = 'coordinate',
                               output_filename = output_file_name,
                               ending_pattern  = ending_pattern)
    basename = sample_name + '.myco'
    _build_BAM_index(basename)
    _extract_samtool_stats(output_file_name)
    _run_BWA_reversed_mapping( output_file_name,
                               suffix           = '.myco.bam',
                               bwa_index_file   = bwa_index_file,
                               output_file_name = None,
                               threads          = THREADS)
    myco_rev_bam_file = sample_name + '.rev.bam'
    print('myco_rev_bam_file ', myco_rev_bam_file)
    _extract_samtool_stats(myco_rev_bam_file, suffix = '.rev.bam')
    
    return None

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
@parameters 
    1. choice(String)
               the choice defined when performing the script.
@return(string) :
    String. from inner dictionary. This will initilize the param
    for the next step.
"""
def switch_genome_build(choice):
    
    return {
        'hg38': HG38_UCSC_GENOME,
        'hg19': HG19_UCSC_GENOME,
        'mm10': MM10_UCSC_GENOME,
        'mm9' : MM9_UCSC_GENOME,
        'rn6' : RN6_UCSC_GENOME
    }[choice]

def switch_BWA_index(choice):
    
    return {
        'hg38': BWA_INDEX_HG38_PATH,
        'hg19': BWA_INDEX_HG19_PATH,
        'mm10': BWA_INDEX_MM10_PATH,
        'mm9' : BWA_INDEX_MM9_PATH,
        'rn6' : BWA_INDEX_RN6_PATH
    }[choice]

def switch_strandness(choice):
    
    return {
        'NONE': 'NONE',
        'FR'  : 'FIRST_READ_TRANSCRIPTION_STRAND',
        'RF'  : 'SECOND_READ_TRANSCRIPTION_STRAND'
    }[choice]

def switch_ref_flat_picard(choice):
    return {
        'hg38'  : REFFLAT_HG38_UCSC_PICARD ,
        'hg19'  : REFFLAT_HG19_UCSC_PICARD ,
        'mm10'  : REFFLAT_MM10_UCSC_PICARD,
        'mm9'   : REFFLAT_MM9_UCSC_PICARD ,
        'rn6'   : REFFLAT_RN6_UCSC_PICARD 
    }[choice]

def switch_ribo_interval(choice):
    return {
        'hg38'  : RIBO_INTERVAL_LIST_HG38_PICARD ,
        'hg19'  : RIBO_INTERVAL_LIST_HG19_PICARD ,
        'mm10'  : RIBO_INTERVAL_LIST_MM10_PICARD,
        'mm9'   : RIBO_INTERVAL_LIST_MM9_PICARD ,
        'rn6'   : RIBO_INTERVAL_LIST_RN6_PICARD
    }[choice]

'''
@aim:
   this will change the working dir using the outside defined param
   which is required and used in the next step. the working dir
   should have the write/read model. Otherwise, will throw out the
   error.

@parameters
    1. working_dir(String): the working dir used in the next.
                            all program the temporary results or final
                            output will be saved in this dir.

@return(None)
@update 2018-04-19

'''


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



'''
@aim
    Generating an MD5 checksum of a file
    https://stackoverflow.com/questions/3431825/
    generating-an-md5-checksum-of-a-file
    
    Write MD5 hashes to file for all files in a directory tree
    https://codereview.stackexchange.com/questions/133859/
    write-md5-hashes-to-file-for-all-files-in-a-directory-tree

    https://www.joelverhagen.com/blog/2011/02/md5-hash-of-file-in-python/
@parameter
    filename(string): the file name with path
@return (string)
    the finger print of a single file
@update  2018-03-19
'''
def _md5sum(filename):

    hash_md5 = hashlib.md5()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

'''
@aim 
    to generate the each fingerprints of the file
    using the function above.
    the fingerprints are md5sum.
@parameter
    list
@return (dict) :
    the key is the filename with path
    the value is its fingerprint
    the fingerprints
@update 2018-03-19
'''
    
def _generate_raw_data_fingerprints(filenames):
    assert isinstance(filenames,list)
    assert len(filenames) != 0
    fingerprints_dict = dict()
    for file in filenames:
        assert os.path.isfile(file) 
        fingerprint = md5sum(file)
        fingerprints_dict[file] = fingerprint
    return fingerprints_dict

'''
@aim 
    this function will check the validity of the raw data.
    if the data will fulfill the naive quality level and 
    whole and accuracy.
@param
    1. prints_dict(dict): finger prints in the dictionary form
@return
    save the finger prints in a local file (txt)
'''

def save_raw_data_fingerprints( filenames,
                                output_file = 'fingerprints.txt'):
    prints_dict = _generate_raw_data_fingerprints(filenames)
    with open(output_file, 'w') as fh:
        for key in prints_dict:
            fh.write('%s\t%s\n' % (key, prints_dict[key]) )
    
    return output_file

def check_raw_data_fingerprints():
    return None

def compute_running_time(start, stop):
    computation_time = stop - start
    print('computation_time is %s in second\n' % computation_time )
    return None

def read_fastQC_report():
    return None

def read_Picard_report():
    return None


'''
@aim:
   this is the top function and will perform the task.

@parameters
    1. params_object: the param object.

@return(None)
@update 2018-03-09

'''
def perform_mRNA_QCtask():
    param_dict = vars(param_parser.parse_args())
    
    #---
    # now get all required parameters from script inputs
    #---
    
    STRANDNESS          = switch_strandness( param_dict['strandness'] )
    BWA_INDEX_PATH      = switch_BWA_index( param_dict['genome_build'] )
    REFERENCE_GENOME    = switch_genome_build( param_dict['genome_build'] )
    REF_FLAT            = switch_ref_flat_picard( param_dict['genome_build'] )
    RIBOSOMAL_INTERVALS = switch_ribo_interval( param_dict['genome_build'] )
    file_suffix         = param_dict['name_pattern']
    QC_type             = param_dict['QC_type']
    data_path           = param_dict['data_path']
    library_model       = param_dict['library_model']
    THREADS             = param_dict['threads']
    working_path        = param_dict['working_path']
    
    whole_data_names    = get_raw_data_names( 
                              data_path, 
                              ending_pattern = file_suffix)
    set_working_path(working_path)
    #run_FASTQC(whole_data_names)
    
    #run_multiThreads_FASTQC( whole_data_names, 
    #                         threads      = THREADS,
    #                         output_dir   = 'fastqc.results')
    
    if library_model == 'PE':
        read1_list, read2_list = split_PairEnd_files(whole_data_names)
        whole_data_names       = read1_list
        for i in range(len(read1_list)):
            sample_name = run_BWA_with_limit_memory( 
                                read1_list[i], 
                                read2_list[i],
                                library_model   = library_model,
                                bwa_index_file  = BWA_INDEX_PATH,
                                threads         = THREADS,
                                output_filename = None,
                                ending_pattern  = file_suffix)
            run_PICARD_QC_modules( 
                    sample_name,
                    ref_genome      = REFERENCE_GENOME,
                    ref_flat        = REF_FLAT,
                    ribo_annotation = RIBOSOMAL_INTERVALS,
                    strandness      = STRANDNESS)
    elif library_model == 'SE':
        for i in range(len(whole_data_names)):
            sample_name = run_BWA_with_limit_memory( 
                                whole_data_names[i], 
                                library_model   = library_model,
                                bwa_index_file  = BWA_INDEX_PATH,
                                threads         = THREADS,
                                output_filename = None,
                                ending_pattern  = file_suffix)
            run_PICARD_QC_modules( 
                    sample_name,
                    ref_genome      = REFERENCE_GENOME,
                    ref_flat        = REF_FLAT,
                    ribo_annotation = RIBOSOMAL_INTERVALS,
                    strandness      = STRANDNESS)
    for raw_data in whole_data_names:
        check_mycoplasma_contamination( 
                    raw_data,
                    read2           = None,
                    library_model   = 'SE',
                    myco_bwa_index_file  = MYCOPLASMA_GENOMES_BWA_INDEX,
                    myco_sam_index_file  = MYCOPLASMA_GENOMES,
                    bwa_index_file       = BWA_INDEX_PATH,
                    sam_index_file       = REFERENCE_GENOME,
                    threads              = THREADS,
                    sorting_method       = 'coordinate',
                    ending_pattern       = file_suffix )
    print('now, we have completed the QC task if no errors are thrown-out!')
    return None

'''
@aim
    to perform the ChIPseq QC procedure
@parameter
@return
@update   2018-03-23
'''
def perform_ChIP_QCtask(params_object):
    return None


"""
@aim   open the debugg model to find and test function
       now the script get the input parameters
       which are specified by the end user.
@param None

@return (list) paired end seq files

cd /home/zhenyisong/data/cardiodata/mycoplasma

rm -rf *.bam *.bai *.txt *.out *fastqc *.stats *.out *.html *.zip

check_mycoplasma_contamination( 'SRR488569.fastq.gz',
                                read2           = None,
                                library_model   = 'SE',
                                myco_bwa_index_file  = MYCOPLASMA_GENOMES_BWA_INDEX,
                                myco_sam_index_file  = MYCOPLASMA_GENOMES,
                                bwa_index_file       = BWA_INDEX_PATH,
                                sam_index_file       = REFERENCE_GENOME,
                                threads         = THREADS,
                                sorting_method  = 'coordinate',
                                middle_name     = '',
                                ending_pattern  = 'fq.gz' )

"""

def debug_model():
    get_test_data_path =  os.path.dirname(sys.argv[0])
    print( '''now you open the debug model, we will 
              use the sample data to simulate the results''')
    param_dict = vars(param_parser.parse_args())
    #---
    # get all required parameter to perform debug model
    #---
    STRANDNESS          = switch_strandness( param_dict['strandness'] )
    BWA_INDEX_PATH      = switch_BWA_index( param_dict['genome_build'] )
    REFERENCE_GENOME    = switch_genome_build( param_dict['genome_build'] )
    REF_FLAT            = switch_ref_flat_picard( param_dict['genome_build'] )
    RIBOSOMAL_INTERVALS = switch_ribo_interval( param_dict['genome_build'] )
    file_suffix         = '.downsample.fq.gz'
    QC_type             = param_dict['QC_type']
    library_model       = param_dict['library_model']
    THREADS             = param_dict['threads']
    working_path        = param_dict['working_path']
    data_path           = None
    if ( STRANDNESS != 'NONE'                  and  
         param_dict['genome_build'] == 'mm10'  and 
         library_model    == 'PE'):

        data_path           = get_test_data_path + '/test/SRP082391'
    elif ( STRANDNESS    == 'NONE'              and 
           param_dict['genome_build'] == 'mm10' and 
           library_model    == 'PE'):
        data_path           = get_test_data_path + '/test/SRP109298'
    elif ( STRANDNESS       != 'NONE'           and 
           param_dict['genome_build'] == 'hg38' and 
           library_model    == 'PE'):
        data_path           = get_test_data_path + '/test/SRP124631'
    elif ( STRANDNESS       != 'NONE'            and 
           param_dict['genome_build'] == 'rn6'   and 
           library_model    == 'PE'):
        data_path           = get_test_data_path + '/test/SRP074376'
    else:
        print('no sample data set for your request!')
        sys.exit(0)

    whole_data_names    = get_raw_data_names( 
                              data_path, 
                              ending_pattern = file_suffix)
    set_working_path(working_path)

    #print(whole_data_names)
    
    if library_model == 'PE':
        read1_list, read2_list = split_PairEnd_files(whole_data_names)
        for i in range(len(read1_list)):
            sample_name = run_BWA_with_limit_memory( 
                                read1_list[i], 
                                read2_list[i],
                                library_model   = library_model,
                                bwa_index_file  = BWA_INDEX_PATH,
                                threads         = THREADS,
                                output_filename = None,
                                ending_pattern  = file_suffix)
            run_PICARD_QC_modules( 
                    sample_name,
                    ref_genome      = REFERENCE_GENOME,
                    ref_flat        = REF_FLAT,
                    ribo_annotation = RIBOSOMAL_INTERVALS,
                    strandness      = STRANDNESS)
    elif library_model == 'SE':
        for i in range(len(whole_data_names)):
            sample_name = run_BWA_with_limit_memory( 
                                whole_data_names[i], 
                                library_model   = library_model,
                                bwa_index_file  = BWA_INDEX_PATH,
                                threads         = THREADS,
                                output_filename = None,
                                ending_pattern  = file_suffix)
            run_PICARD_QC_modules( 
                    sample_name,
                    ref_genome      = REFERENCE_GENOME,
                    ref_flat        = REF_FLAT,
                    ribo_annotation = RIBOSOMAL_INTERVALS,
                    strandness      = STRANDNESS)
    print('now, we have completed the debugging if no errors are thrown-out!')
    return None


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
   
start_time = timeit.default_timer()

if vars(param_parser.parse_args())['debug']:
    debug_model()
else:
    perform_mRNA_QCtask()
    sleep(10)

stop_time  = timeit.default_timer()

compute_running_time(start_time, stop_time)

sys.exit(0)

