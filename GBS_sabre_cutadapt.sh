#!/bin/bash -e
#SBATCH --job-name=GBS_preproc      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=6:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=1      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis   # directory where you run the job
#SBATCH --output=%x-%j.out      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# This script was adapted for use in NeSI with a SLURM system, on 20.11.18
# loading necessary modules:
module load cutadapt/1.16-gimkl-2017a-Python-3.6.3
sabre=/nesi/project/uoo02327/denise/sabre-master/sabre

# defining variables:
datadir=/nesi/nobackup/uoo02327/denise/ModPop_analysis
datafile=SQ0501_S6_L006_R1_001.fastq.gz

# setting up a logfile
echo "Logfile for GBS preprocessing run on "`date` > preprocess.log
logfile=${datadir}/preprocess.log

echo "Demultiplexing with sabre" >> $logfile

mkdir demultiplexed
cd demultiplexed/

# command to run sabre for single end data, uncomment the option needed: -m 1 allows for 1 mismatch in barcode
#sabre se -f ${datadir}/data.fq -b ${datadir}/barcodes.txt -u unknown_barcode.fq > sabre_summary.txt
$sabre se -f ${datadir}/${datafile} -m 1 -b ${datadir}/barcodes.txt -u unknown_barcode.fq > sabre_summary.txt

echo "Demultiplexing done"`date` >> $logfile

# removing unnecessary files for later steps
rm GBSNEG*.fq
rm unknown_barcode.fq

# calling my sample list for the next loop from the demultiplexed files here
samplist=`ls -1 *.fq | sed 's/.fq//'`

cd ..

#####################################

echo "Filtering and trimming with cutadapt" >> $logfile

echo "Preprocessing summary" > preprocessing_summary.txt
summary=${datadir}/preprocessing_summary.txt
echo "Sample"$'\t'"Demultiplexed_reads"$'\t'"Filtered_reads"$'\t'"Trimmed reads" >> $summary

mkdir filtered
mkdir trimmed

for samp in $samplist

do

  echo "Processing $samp "`date` >> $logfile

  cd filtered

  echo -n "$samp"$'\t' >> $summary
  echo -n `grep '^@' ../demultiplexed/${samp}.fq | wc -l`$'\t' >> $summary

  ## this is to select only the reads that begin with the proper enzyme restriction site
  grep -B1 -A2 '^TGCAG' ../demultiplexed/${samp}.fq | sed '/^--$/d' > ${samp}.fq

  echo -n `grep '^@' ${samp}.fq | wc -l`$'\t' >> $summary

  cd ..

  #####

  cd trimmed/

  # command to run cutadapt
  cutadapt -a file:${datadir}/adaptersSE.fa -m 50 -o ${samp}.fq ../filtered/${samp}.fq >> cutadapt_summary.txt

  ## using an adapter.fasta file with the adapter sequences that come with trimmomatic, checked with Illumina, they should be fine, no quality trimming because BWA does that on its own; min length to keep a read is 50

  grep '^@' ${samp}.fq | wc -l >> $summary

  cd ..

  echo "$samp processed" >> $logfile

done

echo "Filtering and trimming done "`date` >> $logfile

######################################
