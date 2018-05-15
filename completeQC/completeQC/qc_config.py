
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

# test data set location

TEST_DATA_PATH   = '/wa/zhenyisong/sourcecode/core.facility/completeQC'

#--- index end
