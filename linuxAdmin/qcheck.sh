#!/bin/bash

# qsub /wa/zhenyisong/sourcecode/core.facility/linuxAdmin/qccheck.sh

# please see more reference at this page
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
# command-line interpreter/shell �� see below for details
#---

#$ -S /bin/bash

#---
# The name of the job. The name should follow the  "name"
# definition  in sge_types(1).  Invalid job names will be
# denied at submit time.
#---

#$ -N chenxi.QC

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


#---
# Mail at beginning/end/on suspension
#---

#$ -m ea

#---
# Send mail to these users
#---

#$ -M zhenyisong@hotmail.com
#$ -wd /home/zhenyisong/data/results/chenlab/upload/Cleandata
#$ -j yes
#$ -o job.log
#$ -e error.log

source activate biotools
unset PYTHONPATH
quality_control_industry.py -n '.fq.gz' -g 'hg38' -l 'PE' -s 'NONE'

source deactivate biotools