#!/bin/bash -e
#SBATCH --job-name=GBS_GATK      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=24:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=8      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis   # directory where you run the job
#SBATCH --output=%x-%j.out      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# This script was adapted for use in NeSI with a SLURM system, on 22.11.18
# loading necessary modules:
module load GATK/3.8-1
gatk=/opt/nesi/mahuika/GATK/3.8-1/GenomeAnalysisTK.jar

# input files and working directory
datadir=/nesi/nobackup/uoo02327/denise/ModPop_analysis
ref=$datadir/pseudochromosomes.fasta
bamdir=$datadir/alignment
logfile=$datadir/gatk_varcaller.log
vcfdir=$datadir/gatk_vcf


cat samplelist.txt | while read name
do
  echo "Calling variants for $name "$(date) >> $logfile
  # this produces variants for each .bam alignment
  java -jar $gatk -T HaplotypeCaller -R $ref -I $bamdir/realigned_$name\.bam --emitRefConfidence GVCF -o $vcfdir/$name\_raw.snps.g.vcf

done

echo "Setting up the .vcf list "$(date) >> $logfile
for name in $(cat samplelist.txt)
do
  # this loop prepares the input list for GenotypeGVCFs (--variant file1.g.vcf --variant file2.g.vcf ... --variant fileN.g.vcf)
  variant=' --variant '$vcfdir/$name\_raw.snps.g.vcf
  #echo $variant
  allVariants="${allVariants}$variant"
done

echo $allVariants
echo $allVariants >> $logfile

echo "Merging variant calls "$(date) >> $logfile
# this is the GATK command that unifies the individual vcf calls:
java -jar $gatk -T GenotypeGVCFs -R $ref -o GATK_output.vcf $allVariants
