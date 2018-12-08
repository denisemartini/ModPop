#!/bin/bash -e
#SBATCH --job-name=GBS_samtools      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=08:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=4      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis/   # directory where you run the job
#SBATCH --output=%x-%j.log      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# loading samtols and bcftools
module load SAMtools
module load BCFtools

datadir=/nesi/nobackup/uoo02327/denise/ModPop_analysis
ref=$datadir/pseudochromosomes.fasta

## modified commands to include other FORMAT fields
#e.g. -t DP,AD,INFO/AD

samtools mpileup -min-MQ 20 --min-BQ 20 --BCF --uncompressed --output-tags DP,AD,INFO/AD,INFO/DPR \
--fasta-ref $ref \
--bam-list bamlist.txt | bcftools call \
--format-fields GQ,GP --variants-only --multiallelic-caller \
--output-type v - --output samtools_output.vcf
