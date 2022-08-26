#!/bin/bash -e
#SBATCH --job-name=GBS_map      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=36:00:00         # Walltime (HH:MM:SS)
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
fqdir=$datadir/trimmed
logfile=$datadir/alignment.log

# to clean up as you go
mkdir mapped
mkdir rg
mkdir dedup
mkdir realigned

echo "Starting run "$(date) > $logfile
## first, check if the reference genome has been indexed:
if [ ! -f $ref\.amb ]
then
  bwa index -a bwtsw $ref
  else
    echo "BWA index found" >> $logfile
fi
if [ ! -f $ref\.fai ]
then
  samtools faidx $ref
  else
    echo "Samtools .fai index found" >> $logfile
fi


#Create a Sequence Dictionary if necessary
if [ ! -e $ref\.dict ]
then
  echo "SequenceDictionary of reference does not exist: creating one with picard"
  java -jar $picard CreateSequenceDictionary \
  R=$ref \
  O=$(echo $ref | cut -f 1 -d '.').dict
  else
    echo "SequenceDictionary found" >> $logfile
fi

# creating a summary file of mapping
echo "Mapping summary" > mapping_summary.txt
summary=${datadir}/mapping_summary.txt
echo "Sample"$'\t'"Mapped_reads"$'\t'"%_mapped"$'\t'"Unmapped_reads" >> $summary

# start processing

cat fq_list.txt | while read fq
do
  echo $fq
  echo "$fq started "$(date) >> $logfile
  name=$(echo $fq | cut -f 1,2 -d '_')
  # alignment
  echo "processing alignment "$(date) >> $logfile
  # bwa mem -t 8 genome.fa reads.fastq | samtools sort -@8 -o output.bam -
  bwa mem -t 8 $ref $fqdir/$fq | samtools sort -@8 -o sorted_$name\.bam -
  samtools index sorted_$name\.bam

  # adding mapping stats to summary
  map=$(samtools view -F4 -c sorted_$name\.bam)
  unmap=$(samtools view -f4 -c sorted_$name\.bam)
  total=$(($map + $unmap))
  perc_mapped=`echo "scale=4;($map/$total)*100" | bc`
  echo -n $name$'\t' >> $summary
  echo -n $map$'\t' >> $summary
  echo -n $perc_mapped$'\t' >> $summary
  echo $unmap >> $summary

  # fix read groups
  echo "fixing read groups "$(date) >> $logfile
  rgpu=$(echo $fq | cut -f 1 -d '.' | cut -f 3,4 -d '_')

  java -jar $picard AddOrReplaceReadGroups \
    I=sorted_$name\.bam \
    O=rg_$name\.bam \
    RGID=D00390:320 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=$rgpu \
    RGSM=$name

  mv sorted*.bam mapped
  mv sorted*.bai mapped

  # fixing everything else the way gatk likes it
  echo "marking duplicates "$(date) >> $logfile
  java -jar $picard MarkDuplicates INPUT=rg_$name\.bam OUTPUT=/dev/stdout METRICS_FILE=$name\_metrics.txt | \
  java -jar $picard ReorderSam I=/dev/stdin O=dedup_$name\.bam R=$ref
  samtools index -b dedup_$name\.bam
  mv rg*.bam rg

  echo "realigning indels "$(date) >> $logfile
  java -jar $gatk -T RealignerTargetCreator -R $ref -I dedup_$name\.bam -o $name\.intervals
  java -jar $gatk -T IndelRealigner -R $ref -I dedup_$name\.bam -targetIntervals $name\.intervals -o realigned_$name\.bam
  samtools index -b realigned_$name\.bam
  mv dedup*.bam dedup
  mv dedup*.bai dedup
  mv realigned*.bam realigned
  mv realigned*.bai realigned
  rm *.intervals

  echo "$fq done "$(date) >> $logfile

done
