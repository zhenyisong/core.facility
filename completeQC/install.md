#!/bin/bash

#---
# @author Yisong Zhen
# @since  2018-05-03
# @update 2018-05-15
#---

#---
# download miniconda
#
# get and the required miniconda package
# from its website:
# https://conda.io/miniconda.html
# https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# you have to choose the recent stable
# version using >=python3.6, not the python 2.x or below
# the QC pipeline was implemented using the python3.x as the 
# engine to perform the NGS data analysis.
# However, the choice of which Miniconda is 
# installed only affects the root environment.
#---

```
miniconda_address='https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh'
wget -c ${miniconda_address} -O miniconda.sh
```

#---
# check fingerprint
#
# if susscefully download the right version
# you can check the package md5 sum to verify
# the validity of the installation package.
#---

```
mdsum miniconda.sh
```


#--- 
# install the miniconda
#
# in the installation process
# you should 'YES' and use its
# default setting to continue
# the installation
#---

```
bash miniconda.sh
```

#---
# configure the minicoda
#
# you have to add these channels in the order
# we use the conda mirror in China at USTC
# which provide much better service than
# tuna at Tsignhua.
# https://mirrors.ustc.edu.cn/help/anaconda.html
#---

```
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

conda create --name qcline

```
#---
# now you should unset the PYTHONPATH if your system already
# install the PYTHON somewhere. This variable will complex the 
# situation. You need to unset it.
#---

```
vi ~/.bash_profile
```

#---
# you will add the line
# unset PYTHONPATH
#---

```
source ~/.bash_profile
source activate qcline

conda install bwa picard plumbum samtools gatk4 fastqc bedtools hisat2
```