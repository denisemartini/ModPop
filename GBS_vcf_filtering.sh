#!/bin/bash
# 11.12.18, Denise Martini
# Script to put together a few filtering steps on the common SNP dataset, after pipeline merging and some dataset exploration.
# VCFtools is the main tool used here and it needs to be in your PATH for this to work.

datadir=/nesi/nobackup/uoo02327/denise/ModPop_analysis
datafile=SQ0501_S6_L006_R1_001.fastq.gz

# setting up a logfile
echo "Logfile for GBS preprocessing run on "`date` > preprocess.log
logfile=${datadir}/preprocess.log
