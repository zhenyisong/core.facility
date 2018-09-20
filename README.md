## 1. Brief Introduction 
This archive will collect all my working scripts when I am in the core facility at 
Guozhong B201 ([sklcvd](http://www.sklcvd.org)) from 2017-12. 

My major task is to develop and implement quality control pipelines to evaluate the NGS data generated from various research groups in Guozhong.

I created several folders which will be the deposited to save the raw code 
when I perform the required task. The tracking code will be disseminated to my 
collaborators for their lab log usage.

current task list (need update)
- [x] implement and maintain quality control pipeline for the NGS studies
- [x] Genome wide association studies (GWAS), DNA-seq
- [x] Epigenentic data analysis, ChIP-seq, ATAC-seq etc.
- [x] transcriptomic analysis including mRNA, lncRNA, miRNA, etc.
- [ ] Matabolomic analysis
- [ ] Single-cell analysis, for example, scRNA-seq.
- [ ] Other X-omics data or high-through-put data


## 2. Core Facility Service Items


NGS data include, but not limited to, RNAseq (mRNA, microRNA & lncRNA), ChIPseq and other epigentics data, for example, ChIRPseq or WGBS.

1. data wrangling
    * NGS data Quality Control(QC), including pre-QC and post-QC. At present, I only provide the complete quality control analysis for mRNA sequencing data at present.
     
    * data pre-processing, including these following steps.
      * reads mapping;
      * reads counting;
      * variant calling;
      * peaking calling. 
      These steps cover most NGS data analysis tasks.
    * Microarray data pre-processing, especially about the Affymatrix data analysis.
       
2. data analysis (solid biostatistics training based on Rosner, Fundamentals of Biostatics)
    * biostatistic design for wet procedures
    [Experimental Design and Data Analysis for Biologists]\
    (https://www.amazon.com/Experimental-Design-Data-Analysis-Biologists/dp/0521009766)
       * A/B testing;
       * simple frequency analysis;
       * Randomized Complete (& Balanced Incomplete) Block Designs (RCBD & BIBD);
       * Latin Squares design;
       * Graeco-Latin Squares;
       * Factorial experiments;
       * and other design strategies. 
    * sample size estimation and power analysis.
       * using [PWR package, R](https://cran.r-project.org/web/packages/pwr/index.html)package;
       * using [PROPER package, R](https://www.ncbi.nlm.nih.gov/pubmed/?term=25273110)
    * statistical inference and confidential calculation using re-sample or bootstrap;
    * using bioconductor packages to analyze the NGS data, 
      including the tools for GSEA, GOenrich or other motif 
      analysis procedure. I also provide the customed design for 
      the network analysis using igraphs.

3. NGS data modeling (solid training based on ISLR & 
        [ESL](https://web.stanford.edu/~hastie/ElemStatLearn/) 
        by Trevor Hastie/Robert Tibshirani)
    * Linear regression model & general linear regression model 
      including logistic or lasso.Regression  
      diagnostics and regression model validation/selection.
    * Statistical Machine learning, including Support vector machine or k-Nearest
      Neighbors (KNN), Decision tree models, Bagging and random forest, boosting.
    * Unsupervised learning strategies, including Partitioning clustering, Hierarchical clustering
      Clustering validation and other advanced clustering procedure. Principle component analysis (PCA)
      Singular value decomposition (SVD) Correspondence analysis, PCA, MCA, FAMD, MFA, HCPC etc. 
       
4. NGS data visualization (the fan of [Dr. Hadley Wickham](http://hadley.nz/) in R application)
    * generate R graph using R basic plotting grammar;
    * generate scientific graph using ggplot2 or grid-based packages;
    * scientific graph layout or plot design;
    * interactive plot design and results for web presentation.
       
5. NGS data pipeline implementation
    * Construction of customized analysis pipeline using R engine;
    * Construction of customized analysis pipeline using Python(2/3) engine.
       
6. SQL query and database management
    *  MySQL database; 
    *  Using SQLite to mange the local data and analysis.

7. Rookie training
    * Basic Linux script programming
    * Software development (R/Python)
       
8. Software maintenance
   * reproducible research including construction of computational environmental and Docker image
   * High Performance Cluster. Usage and script coding.

9. Grant application
   * Scientific paper writing
   * Grant application


Supplementary files 1. [service fee from PKU](http://www.bio.pku.edu.cn/displaynews.php?id=7335)