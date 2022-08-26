#!/bin/bash -e
#SBATCH --job-name=GBS_compare_merge      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=00:30:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=1      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=prepost       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis/vcf_filtering   # directory where you run the job
#SBATCH --output=%x-%j.log      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# loading vcftools and path to included perl scripts locations
module load VCFtools
vcftoolsdir=/opt/nesi/mahuika/VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1/bin/

# to compare all the vcfs between them
${vcftoolsdir}vcf-compare -g *_biall_snps.vcf.gz > vcf-compare5.txt

# remember to choose the order of the input files depending on the order of your trust for the genotype calls
${vcftoolsdir}vcf-isec -f -p vcf-isec5 -n +3 \
platypus_output_biall_snps.vcf.gz stacks_output_biall_snps.vcf.gz \
tassel_output_biall_snps.vcf.gz \
samtools_output_biall_snps.vcf.gz GATK_output_biall_snps.vcf.gz

# then I should be able to just concatenate these files
${vcftoolsdir}vcf-concat \
vcf-isec50_1_2_3_4.vcf.gz vcf-isec50_1_3_4.vcf.gz vcf-isec50_2_3_4.vcf.gz \
vcf-isec50_1_2_3.vcf.gz vcf-isec51_2_3_4.vcf.gz vcf-isec50_1_2_4.vcf.gz \
vcf-isec50_1_2.vcf.gz vcf-isec50_1_3.vcf.gz vcf-isec50_1_4.vcf.gz \
vcf-isec50_2_3.vcf.gz vcf-isec50_2_4.vcf.gz vcf-isec50_3_4.vcf.gz \
vcf-isec51_2_3.vcf.gz vcf-isec51_2_4.vcf.gz vcf-isec51_3_4.vcf.gz \
vcf-isec52_3_4.vcf.gz | bgzip -c > common_snps.vcf.gz

# sorting results
${vcftoolsdir}vcf-sort common_snps.vcf.gz | gzip -c > sorted_common_snps.vcf.gz

# filter out individuals with excessive missingness
# vcftools --gzvcf input_file.vcf.gz --missing-indv --out indv_missingness
# hiatus for now, not necessary for these samples

# filter out snps on missing data
vcftools --gzvcf sorted_common_snps.vcf.gz \
--remove-filtered-all --max-missing 0.1 \
--recode-INFO-all \
--recode --out maxmiss90_common_snps
