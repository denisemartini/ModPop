#!/bin/bash -e
#SBATCH --job-name=GBS_biall_filtering      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=00:10:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=1      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=prepost       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis/vcf_filtering   # directory where you run the job
#SBATCH --output=%x-%j.log      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# loading vcftools
module load VCFtools

# vcftools --vcf input.vcf --min-alleles 2 --max-alleles 2 --remove-indels --remove-filtered-all --recode --out out_biallelic
# repeating above command to loop through my files:
for f in $(ls *.vcf)
do
  vcftools --vcf ${f} \
  --min-alleles 2 --max-alleles 2 \
  --remove-indels \
  --remove-filtered-all --recode \
  --recode-INFO-all \
  --out ${f}_biall_snps
done
