#  GBSpipeline.sh
#  
#
#  Created by Denise Martini on 22/09/17.
#
# just the first few steps of the pipeline for now, testing sabre and cutadapt

# set up datadir and file before beginning
datadir=/data/denise/Kaka_GBS
datafile=SQ0501_S6_L006_R1_001.fastq.gz

# setting up a logfile
start=`date`
echo "Logfile for GBS pipeline run on $start" > logfile.txt
logfile=${datadir}/logfile.txt

echo "Demultiplexing with sabre"
echo "Demultiplexing with sabre" >> $logfile

mkdir demultiplexed
cd demultiplexed/

# command to run sabre for single end data, uncomment the option needed: -m 1 allows for 1 mismatch in barcode
#sabre se -f ${datadir}/data.fq -b ${datadir}/barcodes.txt -u unknown_barcode.fq > sabre_summary.txt
sabre se -f ${datadir}/${datafile} -m 1 -b ${datadir}/barcodes.txt -u unknown_barcode.fq > sabre_summary.txt

echo "Demultiplexing with sabre done"
echo "Demultiplexing done" >> $logfile
echo "Cleaning and defining sample list"

# removing unnecessary files for later steps
rm GBSNEG*.fq
rm unknown_barcode.fq

# calling my sample list for the next loop from the demultiplexed files here
samplist=`ls -1 *.fq | sed 's/.fq//'`

cd ..

#####################################

echo "Filtering and trimming with cutadapt"
echo "Filtering and trimming with cutadapt" >> $logfile

mkdir filtered
echo "Filtering summary" > filtered/filtering_summary.txt
mkdir trimmed
echo "Trimming summary" > trimmed/trimming_summary.txt

for samp in $samplist

do

now=`date`
echo "Processing $samp $now" >> $logfile
echo "Filtering $samp"

cd filtered

echo "$samp" >> filtering_summary.txt
grep '^@' ../demultiplexed/${samp}.fq | wc -l >> filtering_summary.txt

## this is to select only the reads that begin with the proper enzyme restriction site
grep -B1 -A2 '^TGCAG' ../demultiplexed/${samp}.fq | sed '/^--$/d' > ${samp}.fq

grep '^@' ${samp}.fq | wc -l >> filtering_summary.txt

cd ..

#####

echo "Trimming $samp"

cd trimmed/

echo "$samp" >> trimming_summary.txt
grep '^@' ../filtered/${samp}.fq | wc -l >> trimming_summary.txt

# command to run cutadapt
cutadapt -a file:${datadir}/adaptersSE.fa -m 50 -o ${samp}.fq ../filtered/${samp}.fq >> cutadapt_summary.txt

## using an adapter.fasta file with the adapter sequences that come with trimmomatic, checked with Illumina, they should be fine, no quality trimming because BWA does that on its own; min length to keep a read is 50

grep '^@' ${samp}.fq | wc -l >> trimming_summary.txt

cd ..

echo "$samp processed"
echo "$samp processed" >> $logfile

done

end=`date`
echo "Filtering and trimming done"
echo "Filtering and trimming done $end" >> $logfile

######################################