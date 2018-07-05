#!/bin/bash
# @author Yisong Zhen
# @since  2018-06-26
# @update 2018-07-04
#---

# qsub /wa/zhenyisong/sourcecode/core.facility/epigeneticX/cardiac_HiC.sh

#---
# @references
#    1. Lajoie BR, Dekker J, Kaplan N. The Hitchhiker's guide to Hi-C analysis:
#       practical guidelines. Methods. 2015 Jan 15;72:65-75. doi:
#       10.1016/j.ymeth.2014.10.031. Epub 2014 Nov 6. PubMed PMID: 25448293;
#    2.Forcato M, Nicoletti C, Pal K, Livi CM, Ferrari F, Bicciato S. Comparison of
#      computational methods for Hi-C data analysis. Nat Methods. 2017
#      Jul;14(7):679-685. doi: 10.1038/nmeth.4325. Epub 2017 Jun 12. PubMed PMID:
#      28604721; 
#    3. Computational tools for Hi-C data analysis. <Quantitative Biology>
#    4. see Biostar: https://www.biostars.org/p/179295/ <TUTORIALS about Hi-C data analysis>
#----


#---
# @raw data
#    GSE96692
#    SRP102202 : PE/Mouse
# 1: Rosa-Garrido M, Chapski DJ, Schmitt AD, Kimball TH, Karbassi E, Monte E,
#    Balderas E, Pellegrini M, Shih TT, Soehalim E, Liem D, Ping P, Galjart NJ, Ren S,
#    Wang Y, Ren B, Vondriska TM. High-Resolution Mapping of Chromatin Conformation in
#    Cardiac Myocytes Reveals Structural Remodeling of the Epigenome in Heart Failure.
#    Circulation. 2017 Oct 24;136(17):1613-1625. doi:
#    10.1161/CIRCULATIONAHA.117.029430. Epub 2017 Aug 11. PubMed PMID: 28802249;
#    PubMed Central PMCID: PMC5648689.
# please read the attached PDF file and read more detailed protocols
#---


#----
# HPC parameters for Sun Grid
#$ -S /bin/bash
#$ -N colonWGBSseq
#$ -V
#$ -w e
#$ -wd /home/zhenyisong/data/cardiodata/SRP102202
#$ -m ea
#$ -M zhenyisong@gmail.com
#$ -j yes
#$ -o job.log
#$ -e error.log
###$ -l h_vmem=16G
#---


#---
# Please read the supplementary data
# and their working code description
#---
source ~/.bashrc
source ~/.bash_profile
source activate biotools
unset PYTHONPATH


raw_data_path='/home/zhenyisong/data/cardiodata/SRP102202'
cd ${raw_data_path}

decompress_sra_data () {
    file_dir=$1
    data_type=$2
    threads=$3

    if [[ -z $data_type ]]; then
        $data_type='PE'
    fi

    if [ -z "$file_dir" ]; then
        file_dir='./'
    fi

    if [[ -z $threads ]]; then
        $threads=2
    fi

    if [[ ${data_type} == 'PE' ]]; then
        find ${file_dir} -type f -name '*.sra' | \
        xargs -n 1 -P $threads -I{} fastq-dump --split-files --gzip {}
    else
        find ${file_dir} -type f -name '*.sra' | \
        xargs -n 1 -P $threads -I{} fastq-dump --gzip {}
    fi
}

decompress_sra_data ./ PE 3



##mm10_igenome='/wa/zhenyisong/reference/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
##
##mm10_BWA_index='/wa/zhenyisong/reference/Mus_musculus/UCSC/mm10/Sequence/BWAIndex'
##raw_data_path='/home/zhenyisong/data/cardiodata/SRP102202'
##threads=2
##
##if [ ! -d "${WGBS_index}" ]; then
##    mkdir -p ${WGBS_index}
##fi
##
##cd ${WGBS_index}
##
##if [ ! -f "genome.fa" ]
##then
##    ln -s ${hg38_igenome} ./
##fi##