#!/bin/bash

# @author Yisong Zhen
# @since  2018-02-13
# @update 2018-05-04
#---


#---
# perform the apspera download
# chmod 744 /home/zhenyisong/data/sourcecode/core.facility/linuxAdmin/aspera.sh
# nohup /home/zhenyisong/data/sourcecode/core.facility/linuxAdmin/aspera.sh &
#---


#----
# this is the NCBI recommended protocol, but failed
# this need to configure well, otherwise will fail again
# the manual
# https://www.biostars.org/p/111040/
# Tutorial: How to download raw sequence data from GEO/SRA
# prefetch -t ascp --ascp-path \
# '~/.aspera/connect/bin|/etc/aspera|~/.aspera/connect/etc/asperaweb_id_dsa.openssh' \
# SRP032656
#---


# 
# need deletion
# 
# 


ascp='/home/zhenyisong/.aspera/connect/bin/ascp'
openssh='/home/zhenyisong/.aspera/connect/etc/asperaweb_id_dsa.openssh'
ftpadd='ftp-private.ncbi.nlm.nih.gov'
output_dir='/home/zhenyisong/data/cardiodata/'

data='/sra/sra-instant/reads/ByStudy/sra/SRP/SRP055/SRP055990/'

$ascp -T -l 200M -k 3 -i $openssh --host=$ftpadd  \
      --user=anonftp --mode=recv -d ${data} ${output_dir}
