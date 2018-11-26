#!/bin/bash -e
#SBATCH --job-name=GBS_STACKS      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=24:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=10      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis   # directory where you run the job
#SBATCH --output=%x-%j.out      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# This script was adapted for use in NeSI with a SLURM system, on 24.11.18
# loading necessary modules:
module load Stacks/2.0b-gimkl-2017a

# input files and working directory
datadir=/nesi/nobackup/uoo02327/denise/ModPop_analysis
bamdir=$datadir/alignment
outdir=$datadir/stacks

ref_map.pl -T 10 --samples $bamdir --popmap $datadir/population.txt -o $outdir -X "populations:--vcf"
