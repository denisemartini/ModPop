#!/bin/bash
exec 1>vcf_postprocess.log 2>&1
# 11.12.18, Denise Martini
# Script to put together a few filtering steps on the common SNP dataset, after pipeline merging and some dataset exploration.
# VCFtools is the main tool used here and it needs to be in your PATH for this to work.


#set working directory and inputs

datadir=/Volumes/osms-anatomy-dept/Users/D/Denise Martini/Denise/ModPop_analysis/vcf_filtering/
vcffile=maxmiss90_common_snps_withID.vcf
popfile=population.txt
IDfile=high_depth_snp_IDs.txt
indv_file=high_missing_samples.txt

cd $datadir
# setting up a logfile
echo "Logfile for vcf postprocessing run on "`date`
logfile=${datadir}/vcf_postprocess.log

# working variables
basename=$(echo $vcffile | sed 's/\.vcf//')

echo "Filtering specific indvs and loci from $vcffile"
# filter out individuals and snps specified in input files
vcftools --vcf ${vcffile} \
--exclude $IDfile \
--remove $indv_file \
--remove-filtered-all --recode \
--recode-INFO-all \
--out ${basename}

echo "Thinning loci from $vcffile"
# filter out loci that probably came from the same read (i.e. within 80bp of each other)
vcftools --vcf ${basename}.recode.vcf \
--thin 80 \
--remove-filtered-all --recode \
--recode-INFO-all \
--out ${basename}_thinned