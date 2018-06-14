## 1. Brief introduction 
This archive will collect all my working scripts when I am in the core facility at 
Guozhong B201 ([sklcvd](http://www.sklcvd.org)) from 2017-12. 

My major task is to develop and implements quality control pipelines to evaluate the NGS data generated 
from various groups in Guozhong.

I created several folders which will be the deposited to save the raw code 
when I perform the required task. The tracking code will be disseminated to my 
collaborators for their lab log usage.

current task list (need update)
- [x] implement and maintain quality control pipeline for the NGS studies
- [ ] Genome wide association studies (GWAS), DNA-seq
- [ ] Epigenentic data analysis, ChIP-seq, ATAC-seq etc.
- [ ] transcriptomic analysis including mRNA, lncRNA, miRNA, etc.
- [ ] Matabolomic analysis
- [ ] Other X-omics data or high-through-put data


## 2. Core Facility Serivce Items


NGS data inlcude, but not limited to, RNAseq (mRNA, microRNA & lncRNA), ChIPseq and other epigentics data,
for example, ChIRPseq or WGBS.

在以下的服务内容中，我们所指的NGS的数据涵盖，但不局限于，以下数据类型RNASeq（mRNA，lncRNA，microRNA），ChIPseq （或其他epigeneticX的变种，比如ChIRPseq，甲基化WGBS），DNAseq（GWAS，类似）等等，高通量的实验数据。

1. data wrangling
    * NGS data Qaunlity Control(QC), inlucding pre-QC and post-QC. At present, I only provide the complete
      quanlity control analysis for mRNA.
     
    * data preprocessing, inlucing these steps, reads-mapping，reads counting, 
      variant calling & peaking calling. These steps cover most NGS data analysis task.
    * Microarray data prerocessing, epsepcky about the Affimatrix data analuysis.
       
2. data analysis
    1) 实验前期可以进行A/B testing（目前多数的分组分析），以及更为复杂的统计实验设计；实验后期统计假设检验，包括非参数和参数检验，如Student-t 或者是ANOVA等等；
    2) 前期预实验阶段的样本的效能计算（sample size and power law）；
    3) 统计推断下的置信区间计算，采用统计模型或者是数据模拟策略（re-sample）；
    4) 探索性数据分析以及具体统计原理的实现，包括通路富集分析（GSEA），GO分析，网络模型的建立。转录元件的预测（motif prediction）。

3. NGS data modeling
    1) 包括各种线性回归模型（普通的linear regression model，以及广义模型， Logistic/LASSO）采用各种机器学习算法对目前的更大规模的数据进行线性的数据建模和预测；
    2) 包括非线性的预测方法，具体包括具体包括神经网（KNN），支持向量机（SVM），等等传统经典方法；这类工作用于生物标记物的寻找（bio-marker）。更广泛的，是生物数据特征的提取；
    3) 除以上监督学习的策略之外，还可以进行无监督学习，包括各种聚类分析和非线性的聚类分析。比如传统的PCA和其他聚类或降维方法。
       
4. NGS data visualization
    * Genenrate R graph using R basic.使用R basic方式进行统计作图，如Box-plot，Histogram和提琴图等等；
    * Generate Scitific graph using ggplot2 or grid package.
    * 科研作图布局和图形工具的开发（包括Figure的复杂布局）；
    * 基于R环境下的Shiny的网络图形展示和交互（PPT交互和网络数据交互）。
       
5.   NGS data pipeline implementation
    1) 课题组所需的个性化的分析流程，采用R作为引擎进行撰写；
    2) 个性化的分析流程和开发，采用Python作为引擎进行撰写。
       
6. SQL query and database management
    *  MySQL database 
    *  Using SQLite to mange the local data and analysis.

7. Rookie training
    * Basic linux scritp programming
    * Biosofitware develpemtn 
       
8. Software maintenance
    * reproducible research inclduing construcution of computtaonal enveimentand and Dokcer image
    * High Performance Cluster. Usage and script.

9. Grant application
   * Scientific paper writing
   * Grant application


Supplemenray files 1. [service fee from PKU](http://www.bio.pku.edu.cn/displaynews.php?id=7335)