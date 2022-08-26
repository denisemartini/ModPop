## script to run the tassel5-GBS2 pipeline
# setting it up to run from within ModPop_analysis/tassel

# set variables
tassel5="/usr/local/tassel-5-standalone/run_pipeline.pl"
keyfile="/data/denise/ModPop_analysis/keyfile.txt"
enzyme="PstI"
fastq="/data/denise/ModPop_analysis"
DB="./kaka.db"
ref="/data/denise/ModPop_analysis/pseudochromosomes.fasta"


# Step 1: GBSSeqToTagDBPlugin

$tassel5 -Xmx12g -Xms12g -fork1 -GBSSeqToTagDBPlugin -c 3 -e $enzyme -i $fastq -db $DB -k $keyfile -kmerLength 80 -minKmerL 20 -mnQS 10 -endPlugin -runfork1

# Step 2: TagExportToFastqPlugin

$tassel5 -Xmx12g -Xms12g -fork1 -TagExportToFastqPlugin -db $DB -o kaka_tags.fq.gz -c 2 -endPlugin -runfork1

# Step 3:BWA_Alignment

bwa-0.7.12-r1044 aln -n 3 -k 1 -t 10 $ref kaka_tags.fq.gz > kaka_tags.sai
bwa-0.7.12-r1044 samse -n 1 $ref kaka_tags.sai kaka_tags.fq.gz > kaka_tags.sam

# Step 4:SAMToGBSdbPlugin

$tassel5 -Xmx12g -Xms12g -fork1 -SAMToGBSdbPlugin -i kaka_tags.sam -db $DB -aProp 0.0 -aLen 0 -endPlugin -runfork1

# Step 5:DiscoverySNPCallerPluginV2

$tassel5 -Xmx12g -Xms12g -fork1 -DiscoverySNPCallerPluginV2 -db $DB -ref $ref -mnLCov 0.1 -mnMAF 0.01 -endPlugin -runfork1

# Step 6:SNPQualityProfilerPlugin

$tassel5 -Xmx12g -Xms12g -fork1 -SNPQualityProfilerPlugin -db $DB -statFile tassel_stats.txt -endPlugin -runfork1

# Step 7:ProductionSNPCallerPluginV2

$tassel5 -Xmx12g -Xms12g -fork1 -ProductionSNPCallerPluginV2 -db $DB -e $enzyme -i $fastq -k $keyfile -kmerLength 80 -mnQS 10 -o tassel_output.vcf -endPlugin -runfork1
