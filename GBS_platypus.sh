#!/bin/sh

#  GBS_platypus.sh
#  
#
#  Created by Denise Martini on 28/09/17.
#
# to run Platypus with set options

## the following could all be moved to a parameters file
# naming input files
bamlist=`ls *.bam`
refdir=/data/denise
refgen=Kaka_GF3.gapfilled.final.fa
nbcor=10
# these parameters come from the FastGBS pipeline for using Illumina data
minMapQual=20       #min bam quality for that alignment
minBaseQual=20      #min fastq quality for that base
minreads=2          #minimum depth to call a SNP
genIndels=1         #generate indels, 1=yes
# naming output files
logplat=platypus_log.txt
outplat=platypus_output     #only the prefix
##


Platypus.py callVariants --bamFiles="${bamlist}" \
--nCPU="${nbcor}" --minMapQual="${minMapQual}" --minBaseQual="${minBaseQual}" \
--minGoodQualBases=5 --badReadsThreshold=10 \
--rmsmqThreshold=20 --abThreshold=0.01 --maxReadLength=250  --hapScoreThreshold=20 \
--trimAdapter=0 --maxGOF=20 \
--minReads="${minreads}" --genIndels="${genIndels}" --minFlank=5 \
--sbThreshold=0.01 --scThreshold=0.95 --hapScoreThreshold=15 \
--filterDuplicates=0 \
--filterVarsByCoverage=0 --filteredReadsFrac=0.7 --minVarFreq=0.002 \
--mergeClusteredVariants=0 --filterReadsWithUnmappedMates=0 \
--refFile="${refdir}"/"${refgen}" \
--logFileName="${logplat}" \
--output="${outplat}".vcf

