# I used the Dockerfile template here
# http://pythonwise.blogspot.com/2015/04/docker-miniconda-perfect-match.html
# and Install my completeQC python package
#---

FROM miniconda3:latest
MAINTAINER Yisong Zhen <zhenyisong@gmail.com>

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN conda update -y conda


# Python packages from conda
RUN conda install bwa picard plumbum samtools fastqc bedtools hisat2

RUN pip install --no-cache-dir --index-url https://test.pypi.org/simple/ completeQC

