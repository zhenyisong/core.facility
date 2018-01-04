#!/bin/bash

#---
# @author Yisong Zhen
# @since  2018-01-01
# @update 2018-01-03
# SGE setting
#---


# please see more reference at this page
#
# 
# http://gridscheduler.sourceforge.net/htmlman/manuals.html
# http://bioinformatics.mdc-berlin.de/intro2UnixandSGE/
# sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html
# 
# http://talby.rcs.manchester.ac.uk/~ri/
# _linux_and_hpc_lib/sge_intro.html
#---

#---
# These two lines tell SGE to run the script using the BASH 
# command-line interpreter/shell ¡ª see below for details 
#---

#$ -S /bin/bash

#---
# The name of the job. The name should follow the  "name"
# definition  in sge_types(1).  Invalid job names will be
# denied at submit time.
#---

#$ -N Yisong.Aligner

#---
# The -V option to qsub states that the job should 
# have the same environment variables as the 
# shell executing qsub (recommended)
#---

#$ -V

#---
# Specifies a validation level applied to the job  to  be
# submitted  (qsub,  qlogin,  and  qsh)  or the specified
# queued job (qalter). 
#---

#$ -w e

#---
# Execute  the  job  from  the  directory  specified   in
# working_dir.
#---

#$ -wd /home/zhenyisong/data/data/cellChIP/heart

#---
# Mail at beginning/end/on suspension
#---

#$ -m be

#---
# Send mail to these users
#---

#$ -M zhenyisong@gmail.com

#---
# Available for qsub, qsh, qrsh, qlogin and qalter only.
# Specifies whether or not the standard error  stream  of
# the job is merged into the standard output stream.
# If both the -j y and the -e options  are  present,  Sun
# Grid Engine sets but ignores the error-path attribute.
#---

#$ -j yes
#$ -o /wa/zhenyisong/data/log/job.log
#$ -e /wa/zhenyisong/data/log/error.log


###$ -l h_vmem=8G