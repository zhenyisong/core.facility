#!/bin/bash

#---
# @author Yisong Zhen
# @since  2018-12-07
# @update 2018-12-07
# SGE setting
#---



#$ -S /bin/bash

#---
# The name of the job. The name should follow the  "name"
# definition  in sge_types(1).  Invalid job names will be
# denied at submit time.
#---

#$ -N Yisong.HPC

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

#$ -wd /wa/zhenyisong/temp

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
#$ -o /wa/zhenyisong/temp/job.log
#$ -e /wa/zhenyisong/temp/error.log


###$ -l h_vmem=8G
source activate biotools
R CMD BATCH chenxi_chenlab.R
#Rscript a.R
source deactivate biotools