# ModPop

This repo contains the set of scripts and general code used for the processing of GBS data for my modern population study.
All the processing steps (and reasoning behind them) are recorded in the log_ModPop_analysis file.
This file serves as a quick description of the scripts found in the repo.

### ALIGNMENT:

#### GBS_sabre_cutadapt.sh
Script to preprocess a fastq.gz file containing GBS data for all samples.  
Sabre is used to demultiplex the original data, a custom filter is applied to select reads starting with the expected restriction site sequence, then cutadapt is used to trim adapters and short reads. The script is setup to run in NeSI and should be modified to specify the location of working directory, fastq data and sabre and cutadapt executables. Two input files are needed, besides the fastq data: a `barcodes.txt` file for sabre and an `adaptersSE.fa` file for cutadapt; they both need to be in the specified working directory. Output files include a summary file with the number of reads that went through each processing step for each sample, a logfile and the output demultiplexed, filtered and trimmed reads in fastq format.

#### GBS_mapping.sh
To align each preprocessed fastq file to a reference fasta file, assign read groups, remove duplicates, and realign indels. Setup to run in NeSI, the location of working directory, reference fasta file, logfile and executables should be modified to suit. Other parts to check and eventually modify are the read group infos (about halfway through the script). A `fq_list.txt` file need to be present in the working directory, to list the fastq files to process, one per line. This script outputs a clean .bam alignment for each sample, ready to used for variant calling with GATK (the finnicky one) and any other variant caller. The final output files are in the `realigned/` directory. There are intermediate files (.bam alignments) in the `mapped`, `rg` and `dedup` directories. Other outputs include a mapping stats summary, a duplication metrics file for each sample and a run log.

### VARIANT CALLERS:

#### tassel5-GBS2.sh
Script to run all the steps of the Tassel 5 - GBS2 pipeline.  
Setup to run in boros. Should be modified to specify the location of the input files, working directory and Tassel5 executable. The reference fasta file needs to be indexed and named accordingly to the lane and flowcell information present in the keyfile.

#### GBS_GATK.sh
Script to run the GATK variant calling pipeline on NeSI, specifically with GATK version 3.8.1. It will not work with version 4 and above.  
The location of the working directories and inputs should be modified to suit. The reference fasta file needs to be indexed with all the GATK requirements (.dict sequence dictionary, .fai index, etc), which should have been created by the mapping script during realignment. A `samplelist.txt` file needs to be present in the working directory. listing the names of the samples to process, one per line. This scripts outputs individual .vcf files worked by HaplotypeCaller for each sample, in the specified `vcfdir` directory. The final joint calls output, from GenotypeGVCFs, is called GATK_output.vcf and will be in the working directory.
