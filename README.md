## 1. Brief Introduction 
This archive will collect all my working scripts when I am in the core facility at 
Guozhong B201 ([sklcvd](http://www.sklcvd.org)) from 2017-12. 

My major task is to develop and implements quality control pipelines to evaluate the NGS data generated 
from various groups in Guozhong.

I created several folders which will be the deposited to save the raw code 
when I perform the required task. The tracking code will be disseminated to my 
collaborators for their lab log usage.

current task list (need update)
- [x] implement and maintain quality control pipeline for the NGS studies
- [x] Genome wide association studies (GWAS), DNA-seq
- [x] Epigenentic data analysis, ChIP-seq, ATAC-seq etc.
- [x] transcriptomic analysis including mRNA, lncRNA, miRNA, etc.
- [ ] Matabolomic analysis
- [ ] Other X-omics data or high-through-put data


## 2. Core Facility Service Items


NGS data include, but not limited to, RNAseq (mRNA, microRNA & lncRNA), ChIPseq and other epigentics data,
for example, ChIRPseq or WGBS.

1. data wrangling
    * NGS data Quality Control(QC), including pre-QC and post-QC. At present, I only provide the complete
      quality control analysis for mRNA.
     
    * data preprocessing, including these steps, reads-mapping，reads counting, 
      variant calling & peaking calling. These steps cover most NGS data analysis task.
    * Microarray data preprocessing, especially about the Affimatrix data analysis.
       
2. data analysis
    * biostatistic design for the wet
       ** A/B testing
       ** single factor classification (ANOVA)
       ** simple frequency analysis
    * sample size and power law.
    * statistical inference and confidential calculation using re-sample or bootstrap;
    * using bioconductor packages to analyze the data, including the tools, GSEA, GOenrich or
     other motif analysis tools. and network analysis using igraphs.

3. NGS data modeling
    * Linear regression model & general linear regression model including logistic or lasso.
    * Non-linear prediction methods, including SVM or KNN.
    * unsupervised learning strategies. 
       
4. NGS data visualization
    * Generate R graph using R basic.使用R basic方式进行统计作图，如Box-plot，Histogram和提琴图等等；
    * Generate Scientific graph using ggplot2 or grid package.
    * Scientific graph layout or plot
    * interactive plot design and web presentation
       
5.   NGS data pipeline implementation
    * Construction of customized analysis pipeline using R engine.
    * Construction of customized analysis pipeline using Python engine.
       
6. SQL query and database management
    *  MySQL database 
    *  Using SQLite to mange the local data and analysis.

7. Rookie training
    * Basic Linux script programming
    * Software development 
       
8. Software maintenance
   * reproducible research including construction of computational environmental and Docker image
   * High Performance Cluster. Usage and script.

9. Grant application
   * Scientific paper writing
   * Grant application


Supplementary files 1. [service fee from PKU](http://www.bio.pku.edu.cn/displaynews.php?id=7335)