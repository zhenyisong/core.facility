# @author  Yisong Zhen
# @since   2018-10-10
#
# used the Dockerfile template here
# http://pythonwise.blogspot.com/2015/04/docker-miniconda-perfect-match.html
# and Install my completeQC python package


#---
# Run a Docker Image as a Container
#
# https://stackoverflow.com/questions/18497688/run-a-docker-image-as-a-container
#
# first , I should run the docker daemon
# service docker start
# docker build -t completeqc .
# docker run -i -t completeqc quality_control_industry.py -h
##  docker image rm -f completeqc


FROM conda/miniconda3-centos7
MAINTAINER Yisong Zhen <zhenyisong@gmail.com>

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN conda update -y conda


# Python packages from conda
RUN conda install bwa picard plumbum samtools fastqc \
                  bedtools hisat2 psutil setuptools numpy gatk4

RUN pip install --no-cache-dir --index-url https://test.pypi.org/simple/ completeQC
RUN source activate

