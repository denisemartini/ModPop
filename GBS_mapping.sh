#!/bin/bash -e
#SBATCH --job-name=GBS_map      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=12:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=8      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis   # directory where you run the job
#SBATCH --output=%x-%j.out      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# This script was adapted for use in NeSI with a SLURM system, on 21.11.18
# loading necessary modules:
module load BWA
module load SAMtools
module load picard
picard=/opt/nesi/mahuika/picard/2.1.0/picard.jar
module load GATK/3.8-1
gatk=/opt/nesi/mahuika/GATK/3.8-1/GenomeAnalysisTK.jar

# input files and working directory
datadir=/nesi/nobackup/uoo02327/denise/ModPop_analysis
ref=$datadir/pseudochromosomes.fasta
fqdir=$datadir/trimmed/

## first, check if the reference genome has been indexed:
if [ ! -f $ref\.amb ]
	then
	bwa index -a bwtsw $ref
fi
