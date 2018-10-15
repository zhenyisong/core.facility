# @author  Yisong Zhen
# @since   2018-10-15
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
# docker run -i -t completeqc
# quality_control_industry.py -h
##  docker image rm -f completeqc


FROM library/centos:latest
MAINTAINER Yisong Zhen <zhenyisong@gmail.com>

RUN yum install  -y curl
RUN yum install  -y bzip2
# Install miniconda to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN conda update -y conda


# Python packages from conda
RUN conda install bwa picard plumbum samtools fastqc \
                  bedtools hisat2 psutil setuptools numpy gatk4

RUN pip install --no-cache-dir --index-url https://test.pypi.org/simple/ completeQC
#RUN source activate
