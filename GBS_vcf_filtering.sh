#!/bin/bash
exec 1>vcf_postprocess.log 2>&1
# 11.12.18, Denise Martini
# Script to put together a few filtering steps on the common SNP dataset, after pipeline merging and some dataset exploration.
# VCFtools is the main tool used here and it needs to be in your PATH for this to work.

#set inputs

vcffile=maxmiss90_common_snps_withID.vcf
popfile=population.txt
IDfile=high_depth_snp_IDs.txt
indv_file=high_missing_samples.txt

pwd
# setting up a logfile
echo "Logfile for vcf postprocessing run on "$(date)
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

# produce HWE and LD stats for each population in the population.txt file
poplist=$(awk '{print $2}' $popfile | sort -u)

for pop in $poplist
do
  echo "Producing HWE and LD stats for $pop population"
  awk -v var="$pop" '$0~var{print $1}' $popfile > indv.tmp
  vcftools --vcf ${basename}_thinned.recode.vcf \
  --keep indv.tmp \
  --hardy \
  --out ${pop}
  vcftools --vcf ${basename}_thinned.recode.vcf \
  --keep indv.tmp \
  --geno-r2 --min-r2 0.5 --ld-window 20 \
  --out ${pop}
done

echo "Done at "$(date)
