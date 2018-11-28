#!/bin/bash -e
#SBATCH --job-name=GBS_ipyrad      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=04:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=10      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=2              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis   # directory where you run the job
#SBATCH --output=%x-%j.out      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# This script was adapted for use in NeSI with a SLURM system, on 28.11.18
# loading necessary modules:
module load Miniconda3/4.4.10
source activate /nesi/project/uoo02327/programs/miniconda_envs/pyrad

# starting the run, with the demultiplexing and filtering steps
ipyrad -p params-ipyrad.txt -s 12 -c 20 --MPI
# subsetting the run to exclude the negative controls and failed sample
ipyrad -p params-ipyrad.txt -b ipyrad_sub samples.txt
# then running the rest of the pipeline on the remainng samples
ipyrad -p params-ipyrad_sub.txt -s 34567 -c 20 --MPI
