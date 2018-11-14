# @author Yisong Zhen
# @since  2018-11-14
# @references
#    1. Yeo G, Burge CB. Maximum entropy modeling of short sequence motifs with
#    applications to RNA splicing signals. J Comput Biol. 2004;11(2-3):377-94.
#    2. Ito K, Patel PN, Gorham JM, McDonough B, DePalma SR, Adler EE, Lam L, MacRae
#    CA, Mohiuddin SM, Fatkin D, Seidman CE, Seidman JG. Identification of pathogenic
#    gene mutations in LMNA and MYBPC3 that alter RNA splicing. Proc Natl Acad Sci U S
#    A. 2017 Jul 18;114(29):7689-7694.
#
# In response to Jizheng's request, I run the script locally to test it.
# Please read their README file at github : README_regress.score.md
# https://github.com/SplicingVariant/SplicingVariants_Beta
# and the methods section at page 7691-7692 : PMID:28679633
# supporting information at page 1-2        : PMID:28679633

wget https://github.com/SplicingVariant/SplicingVariants_Beta/raw/master/ssfiles.zip
wget https://github.com/SplicingVariant/SplicingVariants_Beta/raw/master/usehg19.tar.gz
wget https://github.com/SplicingVariant/SplicingVariants_Beta/raw/master/usehg38.tar.gz
tar xvzf usehg19.tar.gz & tar xvzf usehg38.tar.gz
unzip ssfiles.zip
cd /wa/zhenyisong/results/songlab/jizheng

__='

R --slave --vanilla --args -mutfile test.mutfile -output mutscore.output \
                           -summarizeresult mutscore \
                           -sjdbout mutsjdbout \
                           -usehg19 usehg19 \
                           -skipRegressScore -skipSRE -refdirectory ssfiles \
                           -skipcDNApos < Regress_Score.v0.97.R
'
R --slave --vanilla --args -mutfile test.mutfile -output mutscore.output \
                           -summarizeresult mutscore \
                           -sjdbout mutsjdbout \
                           -skipRegressScore -skipSRE -refdirectory ssfiles \
                           -skipcDNApos < Regress_Score.v0.97.R