```
#!/bin/bash
```

\@author Yisong Zhen
\@since  2018-05-03
\@update 2018-05-25



- [x] download miniconda first

get and the required miniconda package from [its website](https://conda.io/miniconda.html):
you should choose the recent stable
version default with >=python3.6, not the python 2.x or below
the QC pipeline was implemented using the python3.x as the 
engine to perform the NGS data analysis.
However, the choice of which Miniconda is 
installed only affects the root environment.


```console
miniconda_address='https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh'
wget -c ${miniconda_address} -O miniconda.sh
```


- [x] check fingerprint

If susscefully download the right version
you can check the package md5 sum to verify
the validity of the installation package.

```
mdsum miniconda.sh
```



- [x] install the miniconda

in the installation process
you should 'YES' and use its
default setting to continue
the installation


```
bash miniconda.sh
```


- [x] configure the minicoda

you have to add these channels in the order
we use the conda mirror in China at USTC
which provide much better service than
tuna at Tsignhua. < this is for users in China >
https://mirrors.ustc.edu.cn/help/anaconda.html


```
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

conda create --name <your_env_name>

```

- [x] unset the pythonpath

now you should unset the PYTHONPATH if your system already
install the PYTHON somewhere. This variable will complex the 
situation. You need to unset it. Otherwise, when importing
some libraraies from path lib, there will be errors thrown-out


```
vi ~/.bash_profile
```


you will add the line
unset PYTHONPATH
in your current working enviroment
this way will cancel the conflict when
you have system wide configiration of Python


```
source ~/.bash_profile
source activate <your_env_name>

conda install bwa picard plumbum samtools gatk4 fastqc bedtools hisat2
```

- [x] install the package
now you have to install the my completeQC package 
from pip


```
pip install --no-cache-dir --index-url https://test.pypi.org/simple/ completeQC
```

- [] or uninstall this package

```
pip uninstall completeQC
```


- [x] edit the configiration file

after installation
you should edit the configuration file
which configure the path to several files
including annotaiton files and indexed files
for genome mapping.
the follow conda command will
show the path to the qc_config.py


```
pip show completeQC
```



extract your packge installation path:

```
cd ${Location}/completeQC
```


you need to change the paths in this configuration file.
including the mapping index files for BWA or HISAT2
and other annotaiton files for PICARD.
Following the comments to edit the required file paths.

```
vi qc_config.py
```



