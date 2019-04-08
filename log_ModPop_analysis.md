## Modern Population study - analysis log
#### Denise Martini, 17.11.18
---

All the analysis is done in the ModPop_analysis directory in boros/nesi and moved to HCS for storage

#### What needs to be done:
- [x] Fixing the input files (genome, keyfiles, scripts, etc)
- [x] Rerun variant callers, specifically realignments and then:
  - [x] platypus
  - [x] stacks
  - [x] tassel 5
  - [x] GATK
  - [x] <del>ipyrad</del>
  - [x] samtools
- [x] Filter and merge results, vcftools and VennDiagram
- [ ] Population structure tests:
  - [x] admixture
  - [x] dapc in adegenet
  - [x] <del>tree in treemix</del>
  - [ ] modeling in dadi
  - [x] tree in BEAST
  - [x] migration surfaces with EEMS
- [x] Stats for selection outliers (Fst, Tajima's D, etc)
- [x] Environmental correlations

###### _Technical note_
Most of the first part of the analysis has already been tested before, and scripts/tips are available from the Kaka_GBS directory. Parts yet to test from that part are ipyrad, treemix, dadi.  
One thing that would be worth investigating is if I can run most of this in NeSI this time around, and if it would be faster that way. From a quick check, it should be easy enough to run the alignments, stacks, ipyrad, GATK in NeSI, since they are already installed in there. For the demultiplexing, sabre would need to be added. Platypus and Tassel would also need to be installed, and while Platypus was easy enough (but also fast enough that it would probably be easy to run in boros), Tassel was a nightmare for Hugh the first time around, so I might just run it in boros either way. It was slow, so should be started as quickly as possible. Admixture would probably be an easy install, most of the rest runs in R. Need to investigate dadi.

<!-- TOC -->

- [Modern Population study - analysis log](#modern-population-study---analysis-log)
  - [Denise Martini, 17.11.18](#denise-martini-171118)
  - [What needs to be done:](#what-needs-to-be-done)
  - [Fixing inputs](#fixing-inputs)
  - [Variant Calling](#variant-calling)
    - [TASSEL5](#tassel5)
    - [REALIGNMENT](#realignment)
    - [GATK](#gatk)
    - [PLATYPUS](#platypus)
    - [STACKS](#stacks)
    - [IPYRAD](#ipyrad)
    - [SAMTOOLS](#samtools)
  - [Variant Filtering](#variant-filtering)
  - [Population Structure](#population-structure)
    - [ADMIXTURE](#admixture)
    - [DAPC/ADEGENET](#dapcadegenet)
    - [TREEMIX](#treemix)

<!-- /TOC -->

#### Fixing inputs
###### 19.11.18

First of all, I need to copy all the starting input files I need in the new directory for the analyses.
``` bash
cd /Volumes/osms-anatomy-dept-1/Users/D/Denise\ Martini/Denise/
cp Kaka_assembly/pseudochromosomes.fasta ModPop_analysis
cp Kaka_GBS/SQ0501_S6_L006_R1_001.fastq.gz ModPop_analysis/
cp Kaka_GBS/keyfile.txt ModPop_analysis/
```
Copying everything in boros  
`scp -r ModPop_analysis/ boros:/data/denise/`

Then to modify the names in the fasta file, will need some manipulation
``` bash
cd ModPop_analysis/
sed -i 's/\(>Pseudo\)NC_[0-9]*.[0-9]_/\1_/' pseudochromosomes.fasta
sed -i 's/Pseudo_chromosome/ps_chr/' pseudochromosomes.fasta
sed -i 's/Pseudo_mitochondrion/ps_mito/' pseudochromosomes.fasta
sed -i 's/Superscaffold/un_ssc_/' pseudochromosomes.fasta
sed -i 's/scaffold/un_sc_/' pseudochromosomes.fasta
```

#### Variant Calling
##### TASSEL5
###### 19.11.18

Moving the script to run Tassel5 into the repo and fixing it  
`cp Kaka_GBS/GBS_scripts/tassel5-GBS2.sh ModPop_analysis/ModPop_repo`  
Then, moving to boros and running it  
`scp ModPop_analysis/ModPop_repo/tassel5-GBS2.sh boros:/data/denise/ModPop_analysis/tassel`
``` bash
screen -S TASSEL
cd /data/denise/ModPop_analysis/tassel
bash tassel5-GBS2.sh 2>&1 | tee tassel.log
```
I had forgotten a detail: to run in tassel, the fastq file must be named like the flowcell and lane numbers found in the keyfile, so fixing that before rerunning the above:
```bash
cd ..
mv SQ0501_S6_L006_R1_001.fastq.gz CB67BANXX_6_fastq.gz
cd tassel
bash tassel5-GBS2.sh 2>&1 | tee tassel.log
```
###### 20.11.18

I had forgotten another important detail: Tassel runs BWA at some point for the alignment to the reference .fasta genome, and it assumes that this is indexed already. In this case it wasn't, so it failed to locate the index and did not run the alignment, with obvious consequences for the rest of the pipeline. So:
```bash
cd ..
/usr/local/bwa-0.7.17/bwa index pseudochromosomes.fasta
cd tassel
bash tassel5-GBS2.sh 2>&1 | tee tassel.log
```
###### 21.11.18

Turns out that Tassel5 does weird things with the names of the chromosomes in the reference fasta file, that makes it so it doesn't recognise them anymore later on in the analysis...refer to this link https://groups.google.com/forum/#!topic/tassel/Nd63V-3N-mw  
Anyway, in order to fix this, I will just rename the reference fasta in a way that exclude the word "chr" and that should be enough. Lucky that I started with this variant caller, or I would have had to redo all the others later.
```bash
cd ..
sed -i 's/ps_chr/ps_ch/' pseudochromosomes.fasta
/usr/local/bwa-0.7.17/bwa index pseudochromosomes.fasta
cd tassel
rm *
bash tassel5-GBS2.sh 2>&1 | tee tassel.log
```
_N.B. This time the run finished fine, it took about ~4hrs._

###### 24.11.18

Moving results back to HPC, in the ModPop_analysis directory.
```bash
cd /data/denise/ModPop_analysis/
cp tassel/tassel_output.vcf .
tar -zcvf tassel.gz tassel/
rm -r tassel/
```
`scp boros:/data/denise/ModPop_analysis/tassel.gz .`


##### REALIGNMENT
###### 20.11.18

This will be setup in NeSI, it should be faster there.
```bash
cd /nesi/nobackup/uoo02327/denise
mkdir ModPop_analysis
cd ModPop_analysis/
```
Moving required files to NeSI:  
`scp ModPop_analysis/SQ0501_S6_L006_R1_001.fastq.gz mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`  
`scp ModPop_analysis/pseudochromosomes.fasta mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`  
`scp ModPop_analysis/keyfile.txt mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`

To run the preprocessing (demultiplexing and trimming) I need to install sabre on NeSI:
```bash
cd /nesi/project/uoo02327/denise/
wget https://github.com/najoshi/sabre/archive/master.zip
unzip master.zip
rm master.zip
cd sabre-master/
make
cd ..
```
Installation check, trying to run `./sabre-master/sabre` I get:
```
Usage: sabre <command> [options]

Command:
pe	paired-end barcode de-multiplexing
se	single-end barcode de-multiplexing

--help, display this help and exit
--version, output version information and exit
```
So, the path to sabre to use in my script would be: `/nesi/project/uoo02327/denise/sabre-master/sabre`.  
Now, going back to get a barcodes file ready for sabre, which should be in this format (tab-delimited):
```
CCACCGT   SI-COD01_6_CB67BANXX.fq
GAACAAT   SI-COD02_6_CB67BANXX.fq
AACCGAT   SI-FIO06_6_CB67BANXX.fq
CTGTATG   SI-FIO07_6_CB67BANXX.fq
CCTACAG   SI-FIO08_6_CB67BANXX.fq
```
It should be easy to obtain this from `keyfile.txt`
```bash
cd /nesi/nobackup/uoo02327/denise/ModPop_analysis/
head keyfile.txt
```
```
Flowcell	Lane	Barcode	sample	platename	row	column	libraryprepid	counter	comment	enzyme	species	numberofbarcodes	FullSampleName
CB67BANXX	6	CCACCGT	FT3818	gbs2038	A	1	501			PstI-MspIKaka	96	SI_COD01
CB67BANXX	6	GAACAAT	FT3819	gbs2038	A	2	501			PstI-MspIKaka	96	SI_COD02
CB67BANXX	6	AACCGAT	L34826	gbs2038	A	3	501			PstI-MspIKaka	96	SI_FIO06
CB67BANXX	6	CTGTATG	L40852	gbs2038	A	4	501			PstI-MspIKaka	96	SI_FIO07
CB67BANXX	6	CCTACAG	L42605	gbs2038	A	5	501			PstI-MspIKaka	96	SI_FIO08
CB67BANXX	6	AGCGGTG	L42617	gbs2038	A	6	501			PstI-MspIKaka	96	SI_FIO09
```
I only need fields #3 and #14 for sabre and I can filter them out with `awk`. I also want to skip the header, so I am using `/CB67/` to select the lines that contain that pattern. The `tr` command makes sure that fields are tab-delimited instead of space-delimited. I can then just add the flowcell and lane at the end of the line (`sed` with the `$` symbol), since it is the same for all samples.
```bash
awk -F '\t' '/CB67/{print $3, $14}' keyfile.txt | tr ' ' '\t' >> barcodes.txt
sed -i s'/$/_6_CB67BANXX.fq/' barcodes.txt
```
Moving the script to run sabre and cutadapt into the repo and fixing it for use in NeSI.  
`cp Kaka_GBS/GBS_scripts/GBS_sabre_cutadapt.sh ModPop_analysis/ModPop_repo`  
While fixing the script I realised that I need the adapters file for cutadapt as well, so moving that to the working directory:  
`scp Kaka_GBS/adaptersSE.fa mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`  
Then copying the script over to NeSI:  
`scp ModPop_analysis/ModPop_repo/GBS_sabre_cutadapt.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`  
And running:
```
sbatch GBS_sabre_cutadapt.sh
Submitted batch job 901598
```
The job should be short and not intensive, so the requirements were low, it started right away.  
_N.B. The job ran fine and quite fast, the whole preprocessing took ~2.5hrs. The new summary file also looks good._

###### 21.11.18

I had to fix and index the reference fasta file for Tassel5, so I am transferring that again from boros and into NeSI.  
`scp boros:/data/denise/ModPop_analysis/pseudochromosomes.* ModPop_analysis/`  
`scp ModPop_analysis/pseudochromosomes.* mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`  
Now, for the actual alignment I need to merge together a few scripts from previous tests and adapt them to NeSI.
Since the preprocessing went fine, deleting the unnecessary intermediate files:
```bash
mv demultiplexed/sabre_summary.txt .
mv trimmed/cutadapt_summary.txt .
rm -r demultiplexed
rm -r filtered
```
Creating a new script `GBS_mapping.sh` to align the trimmed sequences to the reference file, add read group information and prepare the alignments for variant calling.
Also preparing a list of the fastq files that need to be processed by the script, in the working directory.
```bash
cd trimmed/
ls * > ../fq_list.txt
cd ..
```
From the above list I am taking off sample `SI_FIO01` because it did not have enough sequences to proceed (pretty much the same as the negative controls, so it is considered a failed sample in the GBS report from AgResearch).
Moving the ready script to NeSI:  
`scp ModPop_analysis/ModPop_repo/GBS_mapping.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`
And running it:
```
sbatch GBS_mapping.sh
Submitted batch job 911985
```
There was error thrown up at the first sample: gatk needs a fasta.fai index for the reference file, so fixing that in the script, adding a separate progress log, and moving it to NeSI again. Also, fixed the sequence dictionary output, it had the wrong name (which GATK would not recognise), and other bits for the log and cleanup.  
`scp ModPop_analysis/ModPop_repo/GBS_mapping.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`
Removing files from the previous run, then restarting it:
```
sbatch GBS_mapping.sh
Submitted batch job 912999
```
If the indel realignment does not take too long, the run should finish all samples within the time limit (36hrs).  
_N.B. The script ran fine, all samples were processed in about 12hrs. The only weird thing I see is in regards to a couple of samples from Zealandia, that seem to have much lower mapping success (in percentage) than the rest. Probably need to check it up?_

###### 22.11.18

The alignment went very well, so I am doing some cleanup and removing intermediate files.
```bash
cat *_metrics.txt >> dedup_metrics_overall.txt
rm *metrics.txt
rm -r mapped/
rm -r rg/
rm -r dedup/
```
Then moving all outputs into a directory named `alignment` and transferring to long storage.  
`scp -r mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/alignment ModPop_analysis`


##### GATK
###### 22.11.18

I have in NeSI everything I need to run the GATK varcaller.  
Creating a new script, called `GBS_GATK.sh`, to run through the individual HaplotypeCaller and the joint genotyping tool of GATK.  
I will need a list of the samples to be processed, in the name only format.
```bash
cut -f1,2 -d '_' alignment/fq_list.txt > samplelist.txt
head -5 samplelist.txt
NI_KAP01
NI_KAP02
NI_KAP03
NI_KAP04
NI_KAP05
```
Moving the new script to NeSI  
`scp ModPop_analysis/ModPop_repo/GBS_GATK.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/`  
Also setting up an output directory `gatk_vcf`, then running  
```
sbatch GBS_GATK.sh
Submitted batch job 929641
```

###### 24.11.18

Unfortunately this time I had not given enough time to the NeSI job (36hrs). The job timed out midway through the individual variant calling of the SI_WES02 sample. To avoid running the whole thing from the beginning, I will modify the script on NeSI, so that it finishes running the missing individuals and then does the joint genotyping.  
Specifically, I am deleting sample SI_WES02 partial output, then creating a samplelist with the samples from this one on:
```bash
rm gatk_vcf/SI_WES02_raw.snps.g.vcf
grep 'SI_WES02' -A 20 samplelist.txt > last_samples.txt
```
Then, on the NeSI script `GBS_GATK.sh` I am switching this new sample list to the previous one, ONLY for the first part of the script.
I am also adding a line to the logfile (to which the rest of the analysis is still going to be appended):  
```bash
echo 'Run interrupted, cleaned and restarted '$(date) >> gatk_varcaller.log
sbatch GBS_GATK.sh
Submitted batch job 970851
```
This time I gave another 24hrs to the job, so hopefully it will get to the end.

###### 26.11.18

The job finished in another 19hrs, so it took overall 55 hrs, fixed the original script accordingly. Specifically, each individual calling took about ~30mins, plus the final joint genotyping took abouth ~12hrs. So cleaning up and moving results to HPC.
```bash
cp GATK/GATK_output.vcf .
tar -zcvf GATK.gz GATK/
rm -r GATK/
```
`scp mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/GATK.gz .`


##### PLATYPUS
###### 24.11.18

To run Platypus in boros (where it has already been installed), I need the alignments from the mapping script that was run in NeSI. Compressing it first, for ease of transfer, from the HPC ModPop_analysis directory.  
```bash
cd alignment
mkdir realigned
mv realigned* realigned/
tar -zcvf ../realigned.gz realigned/
```
`scp realigned.gz boros:/data/denise/ModPop_analysis/`  
In the meantime, adding the platypus script to the repo and fixing it.  
```bash
cp Kaka_GBS/GBS_scripts/GBS_platypus.sh ModPop_analysis/ModPop_repo/
```
Moving to boros, uncompressing alignment directory and running the script:  
`scp ModPop_analysis/ModPop_repo/GBS_platypus.sh boros:/data/denise/ModPop_analysis/`
``` bash
tar -xzvf realigned.gz
screen -S PLATYPUS
mkdir platypus
mv GBS_platypus.sh platypus
cd platypus
bash GBS_platypus.sh 2>&1 | tee platypus_run.log
```

###### 24.11.18

Platypus ran fine on boros, as usual it was very fast and it finished the variant calling in about 1 hour.
Moving results back to HPC, in the ModPop_analysis directory.
```bash
cd /data/denise/ModPop_analysis/
cp platypus/platypus_output.vcf .
tar -zcvf platypus.gz platypus/
rm -r platypus/
```
`scp boros:/data/denise/ModPop_analysis/platypus.gz .`


##### STACKS
###### 24.11.18

I want to setup this in NeSI. Copying a `population.txt` file that I created for previous runs to NeSI.  
`scp Kaka_GBS/popgen/populations/population.txt mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis`  
Wrapping my stacks command in a NeSI script, that I am creating in the repo as GBS_stacks.sh  
```bash
cd ModPop_analysis/ModPop_repo/
nano GBS_stacks.sh
```
Moving the script to NeSI and running it:
```
scp ModPop_analysis/ModPop_repo/GBS_stacks.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis
sbatch GBS_stacks.sh
Submitted batch job 971838
```

###### 26.11.18

The stacks command I used in the past had slightly different options from this version of stacks, so the script failed.  
Fixing that in the GBS_stacks.sh script, then restarting it.  
```
scp ModPop_analysis/ModPop_repo/GBS_stacks.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis
sbatch GBS_stacks.sh
Submitted batch job 1010048
```

###### 27.11.18

The stacks command failed again: I forgot that the name of the alignments must reflect the sample names in the population file.  
So, I decided to quickly fix the names of the alignments for this run:
```bash
cd alignment
for f in $(ls *.ba*); do mv $f $(echo $f | cut -d '_' -f 2-); done
cd ..
```
Also note that stacks wants the output directory to exist already, it does not create it on its own. So `mkdir stacks`.  
And then restart the script:
```
sbatch GBS_stacks.sh
Submitted batch job 1015801
```
_N.B. This script finished really really quickly this time, in about 15 minutes with 10 cores on NeSI. I have fixed the script to reflect that._

###### 28.11.18

Moving results back to HPC, in the ModPop_analysis directory.
```bash
cd /data/denise/ModPop_analysis/
mv population.txt stacks/
cp stacks/populations.snps.vcf ./stacks_output.vcf
tar -zcvf stacks.gz stacks/
rm -r stacks/
```
`scp mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/stacks.gz .`  


##### IPYRAD
###### 28.11.18

I decided a while ago that I wanted to test ipyrad as well, probably in place of samtools, which was not particularly satisfactory the first time around. Problem is that ipyrad is like Tassel5, meaning that it runs the whole pipeline on its own, rather than starting from the alignments. It should be quite fast anyway, and it is already installed on NeSI, but there could be problems of incompatibility at the end, due to the use of different alignments on which the snps are called. Because this problem might arise anyway from Tassel5, I will try using ypirad anyway for now, then decide on what to do once I compare the vcfs.  
I need to fix a parameter file to run ipyrad, but I need to load the conda module and activate the environment in which Hugh installed it on NeSI first. So, to first create a parameters file called params-ipyrad.txt in the ModPop_analysis directory:  
```
module load Miniconda3/4.4.10
source activate /nesi/project/uoo02327/programs/miniconda_envs/pyrad
ipyrad -n ipyrad
```
Inside this parameter file I am going to change all the lines I need for my project:
```
------- ipyrad params file (v.0.7.28)-------------------------------------------
ipyrad                         ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
/nesi/nobackup/uoo02327/denise/ModPop_analysis/ipyrad/   ## [1] [project_dir]: Project dir (made in curdir if not present)
../SQ0501_S6_L006_R1_001.fastq.gz                       ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
../ipyrad_barcodes.txt                               ## [3] [barcodes_path]: Location of barcodes file
                               ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
reference                      ## [5] [assembly_method]: Assembly method (denovo, reference, denovo+reference, denovo-reference)
../pseudochromosomes.fasta                           ## [6] [reference_sequence]: Location of reference sequence file
gbs                            ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TGCAG,                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                          ## [13] [maxdepth]: Max cluster depth within samples
0.85                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
5, 5                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus (R1, R2)
8, 8                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus (R1, R2)
4                              ## [21] [min_samples_locus]: Min # samples per locus for output
20, 20                         ## [22] [max_SNPs_locus]: Max # SNPs per locus (R1, R2)
8, 8                           ## [23] [max_Indels_locus]: Max # of indels per locus (R1, R2)
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus (R1, R2)
0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, v                        ## [27] [output_formats]: Output formats (see docs)
                               ## [28] [pop_assign_file]: Path to population assignment file
```
I will need to modify the barcodes file that I used for the demultiplexing in sabre, specifically to invert the order of columns.
```bash
cp alignment/barcodes.txt ./barcodes.txt
cat barcodes.txt | awk -F '\t' '{print $2, $1}' | tr ' ' '\t' > ipyrad_barcodes.txt
sed -i 's/_6_CB67BANXX.fq//' ipyrad_barcodes.txt
```
Then, creating a list of samples to bring to the end of the analysis, excluding the failed sample (SI_FIO01) and the two negative controls.
```bash
cat ipyrad_barcodes.txt | awk '{print $1}' > samples.txt
nano samples.txt
```
Finally, fixing a script called `GBS_ipyrad.sh`.  
Moving script to NeSI and running it:
```
scp ModPop_analysis/ModPop_repo/GBS_ipyrad.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis
sbatch GBS_ipyrad.sh
Submitted batch job 1051305
```

###### 29.11.18

The ipyrad run unfortunately timed out. It looks like the script is working fine anyway, it took ~1hr to run the first part, that is demultiplexing and filtering on the whole dataset, then the subsetting went fine and the second branch also was running fine.
Because it was not clear at what stage the run had gotten so far, I deleted all the outputs related to the ipyrad_sub part of the script. Then I fixed the script to give it more time (8hrs) and to skip the first two steps (I just commented out the first ipyrad command). And rerunning it:
```
sbatch GBS_ipyrad.sh
Submitted batch job 1059729
```
###### 05.12.18
In the end ipyrad took ~12hrs to comlete the run. There is something confusing in the way the output vcf is formatted: there is first a list of loci that in place of the chrom and position columns are named `locus_####  #` (where # are numbers). Only after these the actual loci with position info start.

I have been examining the ipyrad results and I have come to a few conclusions:
- the reference method is bugged, because the loci that don't align to reference are still assembled denovo and the results are present in the output (when they should not be), they are the records called `locus_###`. So these need to be filtered out from the output .vcf file.
- ipyrad calls "locus" the aligned tag, not the snps...so the output contains ANY variant base from a passing filter LOCUS: as long as there is even only ONE variant in that locus that passes filters, all the other variants (even very crappy ones) ARE present in the output. This is not necessarily a tragedy in itself, because all the crappy ones can be filtered out afterwards. The problem is when bad genotypes are called together at sites where actual alleles are present.
- to add to the above, ipyrad calls as variant sites things that are clearly alignment errors, usually localised at the end of reads, where the trimming has done some damage and left reads with different lengths. These variants are unfortunately not always crappy, meaning that sometimes they have decent depth and would be not easy to filter out. Most of them are tri or tetra-allelic, which is quite uncommon and that's why they are evident from the files, so they would theoretically be filtered out when choosing biallelic sites only.

In an attempt to see if I could clean up some of the bad genotypes at good sites, I used vcftools to filter sites at which an alternate allele is called with a depth less than n reads:
```
vcftools --vcf ipyrad_sub.vcf --non-ref-ac 2 --remove-filtered-all --recode --out ipyrad_output
Outputting VCF file...
After filtering, kept 729794 out of a possible 821695 Sites
vcftools --vcf ipyrad_sub.vcf --non-ref-ac 3 --remove-filtered-all --recode --out ipyrad_output3
Outputting VCF file...
After filtering, kept 286358 out of a possible 821695 Sites
```
###### 06.12.18
On further examination though, this does not seem to eliminate the problem of tri and tetra allelic sites.
```
vcftools --vcf ipyrad_sub.vcf --min-alleles 2 --max-alleles 2
After filtering, kept 473012 out of a possible 821695 Sites
vcftools --vcf ipyrad_output3.recode.vcf --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --out ipyrad_output3_biallelic
After filtering, kept 154786 out of a possible 286358 Sites
```
And even after taking away these sites I see things that I am unhappy with...they might get filtered out later, in the comparison with other programs anyway. I think other than applying all these post filters to this program before comparing it to the others, there is one thing that I could fix from the program parameters themselves, and which is already a requirement in both stacks and tassel pipelines: that I cannot have too many SNPs on the same locus. The default parameter for this is a max of 20 variants per locus and that is really high. So, I am going to try and reduce this parameter to 5 (probably already higher than Tassel, which I think only allows one). I should only need to branch the analysis and run the very last step with this different parameter.
Preparing the new branch:
```
module load Miniconda3/4.4.10
source activate /nesi/project/uoo02327/programs/miniconda_envs/pyrad
ipyrad -p params-ipyrad_sub.txt -b ipyrad_filter

loading Assembly: ipyrad_sub
from saved path: /nesi/nobackup/uoo02327/denise/ModPop_analysis/ipyrad/ipyrad_sub.json
creating a new branch called 'ipyrad_filter' with 93 Samples
writing new params file to params-ipyrad_filter.txt
```
I am directly modifying the parameter (number 22) in the newly created params-ipyrad_filter.txt, bringing it to 5.
I also modified parameter 24, bringing it to 0.8...because it sounds a bit dodgy that a heterozygous sites should not be shared by too many samples: doesn't it depend on the reference? I mean, what if the reference represents only one population, where the homozigosity at that locus is a private allele?
Then I just modified (temporarily for now) the GBS_ipyrad.sh script in NeSI to only do this command:
```
ipyrad -p params-ipyrad_filter.txt -s 7	-c 20 --MPI
sbatch GBS_ipyrad.sh
Submitted batch job 1182413
```
I also reduced the time and memory requirements of the job, since it only needs to do the last step.

###### 07.12.18

The job was still not starting, so I investigated a bit more the NeSI partitions and queuing system and I had forgotten about the prepost partition for small pre and post processing jobs. So, I modified the job requirements in the script to run there:
```
#SBATCH --time=02:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=1      # number of cores per task
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=prepost	  # specify a partition
```
The memory associated with each core is much larger on the prepost partition, so 1 core should be enough. *Edit, one core was not enough, increased to 4 and to 3hrs time*
```
scancel 1182413
sbatch GBS_ipyrad.sh
Submitted batch job 1193080
```
Since this was incurably stuck (NeSI is not happy with me lately, it seems), I ended up transferring everything to boros.
```
scp -r mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/ipyrad/ .
scp mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/params* .
scp mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/ipyrad_barcodes.txt .
```
Then, to actually make it work there I had to fix all the paths in all instances where they appeared, both in the params files and in the assembly (.json) files:
```bash
sed -i 's/nesi\/nobackup\/uoo02327/data/' ipyrad/ipyrad.json
sed -i 's/scale_wlg_nobackup\/filesets\/nobackup\/uoo02327/data/' ipyrad/ipyrad.json
sed -i 's/nesi\/nobackup\/uoo02327/data/' ipyrad/ipyrad_sub.json
sed -i 's/scale_wlg_nobackup\/filesets\/nobackup\/uoo02327/data/' ipyrad/ipyrad_sub.json
sed -i 's/nesi\/nobackup\/uoo02327/data/' params-ipyrad_filter.txt
sed -i 's/nesi\/nobackup\/uoo02327/data/' params-ipyrad_sub.txt
sed -i 's/nesi\/nobackup\/uoo02327/data/' params-ipyrad.txt
```
Finally, I was able to run:  
`module load ipyrad`
`ipyrad -p params-ipyrad_filter.txt -s 7 -c 20 --MPI`
To be honest this was quite fast in boros now (~10 mins), since I am using full 20 cpus because no one else is using it...I need to remember that when NeSI clogs up for me this is still a possibility.

Moving results back to HPC, in the ModPop_analysis directory.
```bash
cd /data/denise/ModPop_analysis/
cp ipyrad/ipyrad_filter_outfiles/ipyrad_filter.vcf ./ipyrad_output.vcf
mv ipyrad_barcodes.txt ipyrad/
mv ipyrad_log.txt ipyrad/
mv GBS_ipyrad.sh ipyrad/
mv params* ipyrad/
tar -zcvf ipyrad.gz ipyrad/
rm -r ipyrad/
```
`scp boros:/data/denise/ModPop_analysis/ipyrad.gz .`  

##### SAMTOOLS
###### 08.12.18
Putting the previously used samtools/bcftools commands in a script, called GBS_samtools.sh and added to repo. Fixing the script for NeSI, adjusting a few quality and output values. Transferring to NeSI, creating a bamlist.txt and starting script.
```
scp ModPop_repo/GBS_samtools.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis
ls alignment/*.bam > bamlist.txt
sbatch GBS_samtools.sh
Submitted batch job 1205275
```
###### 09.12.18
That went quite fast, but the result is a bit overabundant, as usual with samtools, I need to apply a couple of extra filter to the script. Specifically, I will add a filter command at the end of the pipe, to filter based on a minimum depth of 6 reads. It is a bare minimum filter, only to make the pipelines more comparable to start with. I also modified the job requirements, since I verified that the script takes little time.
```
scp ModPop_repo/GBS_samtools.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis
rm samtools_output.vcf
sbatch GBS_samtools.sh
Submitted batch job 1216345
```
Moving results back to HPC, in the ModPop_analysis directory.
```bash
cd /data/denise/ModPop_analysis/vcf_filtering
scp mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/samtools_output.vcf .
```

#### Variant Filtering
###### 06.12.18
I have decided on a few filtering steps that I will be performing on this dataset. First I need to make all datasets a bit more comparable before I compare them. Stacks excludes indels and non-biallelic SNPs by default, so it would probably make sense to do the same on the other pipelines' outputs. This would also already reduce ipyrad's noise, even though there are a lot of multisite loci in that pipeline that would not be in the others, like Tassel...but that's one of the reasons why I am comparing the pipelines in the firstplace. But while I don't need indels and non-biallelic sites for later analysis, some of the multisite loci might be fine, so I am keeping them for now. Platypus and GATK I believe allow some of that, so I will too.
Creating a directory for this in NeSI and HPC and moving all the output files there.
```bash
mkdir vcf_filtering
cd vcf_filtering/
mv ../stacks_output.vcf .
mv ../GATK_output.vcf .
mv ../ipyrad_output.vcf
mv ../platypus_output.vcf .
mv ../tassel_output.vcf .
scp ./* mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/vcf_filtering
```
###### 07.12.18
I also want to count the number of snps called by each pipeline, pre-filtering.
```bash
# first, to remove non-reference loci from ipyrad output
grep -v "locus_" ipyrad_output.vcf > fixed_ipyrad_output.vcf
mv fixed_ipyrad_output.vcf ipyrad_output.vcf
# then, the actual counts
for f in $(ls *.vcf)
do
  echo $f
  grep -v '#' $f | wc -l
done
GATK_output.vcf
71177
ipyrad_output.vcf
178684
platypus_output.vcf
393609
stacks_output.vcf
279162
tassel_output.vcf
120338
```
Putting the biallelic/indels filtering commands in a quick script, to keep the log in the same place for all samples as well.
The script is called `GBS_biall_filtering.sh`. Prepared it to loop through the vcfs in the filtering directory. It is also setup to run in the prepost partition.
Moving to NeSI and starting it.
```bash
scp ModPop_repo/GBS_biall_filtering.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/
sbatch GBS_biall_filtering.sh
Submitted batch job 1196582
```
As expected, the program that lost more loci was platypus.
Quickly fixing the file names:
```bash
for f in $(ls *.recode.vcf)
do
  mv $f $(echo $f | sed 's/.vcf_biall_snps.recode./_biall_snps./')
done
```
And I was forgetting that I will need to do some extra changes to the files before merging, so might as well do them now.
```bash
# to fix reference names in tassel output
sed -i 's/^PS_CH/ps_ch/' tassel_output_biall_snps.vcf
sed -i 's/^UN_SSC/un_ssc/' tassel_output_biall_snps.vcf
# to remove extra samples from tassel output:
module load VCFtools
vcftools --vcf tassel_output_biall_snps.vcf \
--remove-indv "GBSNEG1" --remove-indv "GBSNEG2" --remove-indv "SI_FIO01" \
--remove-filtered-all \
--recode --out tassel_output_biall_snps
mv tassel_output_biall_snps.recode.vcf tassel_output_biall_snps.vcf
rm tassel_output_biall_snps.log
# to reorder sample names in tassel output
for f in $(ls *_biall_snps.vcf)
do
  bgzip $f
  tabix -p vcf $f.gz
done
/opt/nesi/mahuika/VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1/bin/vcf-shuffle-cols -t platypus_output_biall_snps.vcf.gz tassel_output_biall_snps.vcf.gz > fixed_tassel_output.vcf
mv fixed_tassel_output.vcf tassel_output_biall_snps.vcf
# the stacks output also needs to be sorted before tabix works:
/opt/nesi/mahuika/VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1/bin/vcf-sort stacks_output_biall_snps.vcf.gz > sorted_stacks_output.vcf
mv sorted_stacks_output.vcf stacks_output_biall_snps.vcf
# re bgzip and tabix the fixed tassel and stacks files
for f in $(ls *_biall_snps.vcf)
do
  bgzip $f
  tabix -p vcf $f.gz
done
module unload VCFtools
```
###### 08.12.18
Putting together the comparing and merging pipelines code in another script, that I am calling `GBS_pipeline_merge.sh`.
But first, I need to take a look at a few comparisons and see what kind of incompatibilies I have, to decide from which pipelines to pick the common snps.
```bash
module load VCFtools
vcftoolsdir=/opt/nesi/mahuika/VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1/bin/
${vcftoolsdir}vcf-compare -g *_biall_snps.vcf.gz > vcf-compare5.txt
module unload VCFtools
```
There is a problem, when I went to look at the vcf-compare results: ipyrad's call are 98% shared with anyone else. This calls for some investigation, because it sounds very unlikely. My first guess would be that there is something wrong in the way the position is recorded in the vcf output. It looks very much like the positions recorded in ipyrad are somehow shifted from the ones called in all the other programs. The shift is not always the same though, so it is not an easy fix. One thing to figure out is if the different position comes from a difference in the vcf output or in the initial assembly, that seems to be done with a different assembler in ipyrad than in anything else. I went and checked the assembly files that ipyrad outputs at step 3, opened them in tablet and there is nothing wrong with them: the positions where variants are visible correspond exactly with the snp positions called by other programs, like platypus, but ipyrad's call are all wrong...to the point that even the reference base is called wrong. And inconsistently, not like there is a shift that is always the same. I see snps for a sample called at positions where that samples should have a certain depth, following the vcf file, but there are no reads there in the bam files. This is all very suspicious, and I would like to let them know, because they probably have a bug...but not today, I already lost enough time for now. BACK TO VARIANT CALLING, using samtools instead.

###### 09.12.18
Resuming this after going back and working on using samtools instead of ipyrad. So, some clean up:
```bash
cd vcf_filtering
mv ../samtools_output.vcf .
rm ipyrad*
grep -v '#' samtools_output.vcf | wc -l
412080
```
After comparing the vcf files from last time I noticed that recoding with vcftools strips off the info fields from files, unless I specify that it shouldn't do it with the `--recode-INFO-all` option. So, deleting the old filtered files and repeating that again, plus all the following mods on the tassel and stacks files, same steps as above.
```bash
rm *_biall_snps*
scp ModPop_repo/GBS_biall_filtering.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/
sbatch GBS_biall_filtering.sh
Submitted batch job 1216573
```
Quickly fixing the file names:
```bash
for f in $(ls *.recode.vcf)
do
  mv $f $(echo $f | sed 's/.vcf_biall_snps.recode./_biall_snps./')
done
```
And I was forgetting that I will need to do some extra changes to the files before merging, so might as well do them now.
```bash
# to fix reference names in tassel output
sed -i 's/^PS_CH/ps_ch/' tassel_output_biall_snps.vcf
sed -i 's/^UN_SSC/un_ssc/' tassel_output_biall_snps.vcf
# to remove extra samples from tassel output:
module load VCFtools
vcftools --vcf tassel_output_biall_snps.vcf \
--remove-indv "GBSNEG1" --remove-indv "GBSNEG2" --remove-indv "SI_FIO01" \
--remove-filtered-all \
--recode --recode-INFO-all \
--out tassel_output_biall_snps
mv tassel_output_biall_snps.recode.vcf tassel_output_biall_snps.vcf
rm tassel_output_biall_snps.log
# to reorder sample names in tassel output
for f in $(ls *_biall_snps.vcf)
do
  bgzip $f
  tabix -p vcf $f.gz
done
/opt/nesi/mahuika/VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1/bin/vcf-shuffle-cols -t platypus_output_biall_snps.vcf.gz tassel_output_biall_snps.vcf.gz > fixed_tassel_output.vcf
mv fixed_tassel_output.vcf tassel_output_biall_snps.vcf
# the stacks output also needs to be sorted before tabix works:
/opt/nesi/mahuika/VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1/bin/vcf-sort stacks_output_biall_snps.vcf.gz > sorted_stacks_output.vcf
mv sorted_stacks_output.vcf stacks_output_biall_snps.vcf
# re bgzip and tabix the fixed tassel and stacks files
for f in $(ls *_biall_snps.vcf)
do
  bgzip $f
  tabix -p vcf $f.gz
done
module unload VCFtools
```
Now, to quickly check the comparisons:
```bash
module load VCFtools
vcftoolsdir=/opt/nesi/mahuika/VCFtools/0.1.14-gimkl-2017a-Perl-5.24.1/bin/
${vcftoolsdir}vcf-compare -g *_biall_snps.vcf.gz > vcf-compare5.txt
module unload VCFtools
```
I want to quickly check the levels of missingness of the datasets:
```bash
module load VCFtools
for f in $(ls *_biall_snps.vcf.gz)
do
  echo $f
  vcftools --gzvcf $f \
  --max-missing 0.1
done
module unload VCFtools
```
I also ran a few tests to check for the level of missing data in single individuals and that varies a lot between programs, so I will not add it to the merging script for now.
```bash
vcftools --gzvcf stacks_output_biall_snps.vcf.gz --missing-indv --out stac_indv_test
less stac_indv_test.imiss
```
There is really only one problematic individual and I will deal with it later. The merge pipeline script os ready, so:
```bash
scp ModPop_repo/GBS_pipeline_merge.sh mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/
sbatch GBS_pipeline_merge.sh
Submitted batch job 1216736
```
###### 10.12.18
The script ran fine, the final output is a file named `maxmiss90_common_snps.recode.vcf`. Taking a look at this and at the pre-filtering results:
```bash
grep -v '#' maxmiss90_common_snps.recode.vcf | wc -l
161431
zcat common_snps.vcf.gz | grep -v '#' | wc -l
179054
vcftools --vcf maxmiss90_common_snps.recode.vcf --missing-indv --out maxmiss90_indv
```
I have a lot more snps this time compared to last time. The individual missingness looks good, there is only one individual (SI_FIO04) which has ~75% missing data and could be excluded by the analysis. I want to plot the missing data with vcfR, but before that I need to conform the platypus calls in the file to the others, specifically I need a DP field in the FORMAT column. The information is there, but it is called NR. To modify that (and only that):
```bash
cp maxmiss90_common_snps.recode.vcf maxmiss90_common_snps_fixed.vcf
sed -i 's/GQ:NR:NV/GQ:DP:NV/' maxmiss90_common_snps_fixed.vcf
```
All the other pipelines had that information already encoded, so no need to fix anything else.
The other thing I want to do is a VennDiagram plot of the pipeline comparison. I should be able to get all the necessary numbers from the `vcf-compare5.txt` file. I would like to put this in a proper piece of code rather than just copy the numbers over, even if that would be faster right now...  
I am preparing two R scripts:
- one to plot the missingness of the data, in a heatmap and with information per individual (depth and missing data), using vcfR [`Missing_data_vcfR.Rmd`]
- the other one to plot the intersection between pipelines, with VennDiagram and UpSetR [`Intersection_pipelines.Rmd`]
This part of the analysis is in separate reports, in .Rmd format, because I wrote it directly in RStudio.

###### 11.12.18
After the exploratory analysis from yesterday I decided on a couple extra quick filters to apply before proceeding with the HWE and LD tests. Specifically, I noticed a few SNPs with very high depth per sample, which I reckon might be alignments in repetitive regions, a bit dubious.
I have looked into these sites in R, extracting them from the overall dataset and then checking what IDs they have:
```R
library(vcfR)
library(reshape2)
setwd("../vcf_filtering")
vcf <- read.vcfR('maxmiss90_common_snps_fixed.vcf')
dp <- extract.gt(vcf, element = "DP", as.numeric=TRUE)
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
dpf[(dpf$Depth>=1000),] -> dpf1000
unique(dpf1000$Index) -> indexes
indexes
ps_ch_5_71257870 ps_ch_5_71257924 ps_ch_5_71257944 un_sc_4773_6038  un_sc_4773_6011  ps_ch_23_566199  ps_ch_1_9316184  ps_ch_5_3264366 ps_ch_5_3264399
```
It looks like it really is just a handful of loci (close together, so on the same even fewer reads).
Excluding them will be quick. Even if I repeat the search setting the limit at 500 reads per sample I get only ~25 loci, for ~16 reads. Over 100 reads per sample we have 184 loci, again probably from not so many reads. I am honestly inclined to exclude all of these. To do so, though, I am reminded from the vcf file that I will need identifiers for these SNPS, if I want to ask vcftools to exclude them. The way that vcfR gives them IDs automatically is quite neat: they are assigned an ID that is made of `scaffold_position`. That would make it quite easy for me to recognise them even later on. So I think I need to apply another quick fix to this vcf file after all:
```bash
grep "#" maxmiss90_common_snps_fixed.vcf > header.vcf
grep -v "#" maxmiss90_common_snps_fixed.vcf | awk '{OFS="\t"}{$3 = $1"_"$2; print}' > data.vcf
cat header.vcf data.vcf >> maxmiss90_common_snps_withID.vcf
```
It is not going to be that quick on my desktop, but well...
In R I am reading in the file again, with the IDs all set this time (I was having problems with some loci that had a stacks ID before, instead of the position), then I am extracting the IDs of loci with >100 reads per sample and I am writing them in a separate file, one ID per line.
```R
dpf[(dpf$Depth>=100),] -> dpf100
unique(dpf100$Index) -> indexes100
indexes100 <- as.character(indexes100)
lapply(indexes100, write, "high_depth_snp_IDs.txt", append=TRUE, ncolumns=1)
```
Then, putting the rest of the filtering in a `GBS_vcf_filtering.sh` script.  
Running it on my desktop because vcftools is generally quite fast. To test it first:
```bash
vcftools --vcf maxmiss90_common_snps_withID.vcf --chr ps_ch_4A --remove-filtered-all --recode --recode-INFO-all --out test
bash ../ModPop_repo/GBS_vcf_filtering.sh
```
There were a couple of fixes to make (wrong directory setup, error in sample names in population.txt), but it seemed to have worked fine. Set up the correct input file in the script and running again as above.

###### 12.12.18
I am left with 103,563 loci after the filtering, the bulk of it was in the thinning. For the next bit I am joining the per population HWE and LD stats in R, in order to identify the loci that are in disequilibrium in 3 or more populations and eventually filter them out, too. The report of this analysis is in a separate markdown document, called `HWE_LD_analysis.Rmd`.

In R I finished investigating the loci out of HWE and in LD, the results are the files named: `loci_outofHWE.txt`, `loci_inLD.txt` and `loci_inLD_kept.txt`. I am now goin to filter out the selected loci in two rounds from the final vcf file (in case I want a file without HWE but with LD loci later on).
```bash
vcftools --vcf maxmiss90_common_snps_withID_thinned.recode.vcf \
--exclude loci_outofHWE.txt \
--remove-filtered-all --recode \
--recode-INFO-all \
--out maxmiss90_common_snps_HWE
After filtering, kept 92 out of 92 Individuals
Outputting VCF file...
After filtering, kept 103326 out of a possible 103563 Sites
Run Time = 119.00 seconds

vcftools --vcf maxmiss90_common_snps_HWE.recode.vcf \
--exclude loci_inLD.txt \
--remove-filtered-all --recode \
--recode-INFO-all \
--out maxmiss90_common_snps_HWE_LD
```
Finally done with the filtering, too.

#### Population Structure
##### ADMIXTURE
###### 12.12.18
I should be able to use in admixture some .ped plink files, that I can obtain from VCFtools. I am going to set up and run this in boros, because it should be faster and admixture is there already.
```bash
mkdir pop_structure
cd pop_structure/
scp maxmiss90_common_snps_HWE_LD.recode.vcf boros:/data/denise/ModPop_analysis/pop_structure
vcftools --vcf maxmiss90_common_snps_HWE_LD.recode.vcf \
--plink --out maxmiss90_common_snps
```
I might be seeing why in the past I used stacks to convert the files instead: vcftools seems to be doing weird things with my chromosome names. Not sure if that's an issue though. The populations are a problem though, so fixing that:
```bash
sed -i 's/^NI_KAP[0-9]*/KapitiIsland/' maxmiss90_common_snps.ped
sed -i 's/^NI_LBI[0-9]*/LittleBarrierIsland/' maxmiss90_common_snps.ped
sed -i 's/^NI_PUR[0-9]*/Pureora/' maxmiss90_common_snps.ped
sed -i 's/^NI_ZEA[0-9]*/Zealandia/' maxmiss90_common_snps.ped
sed -i 's/^SI_COD[0-9]*/CodfishIsland/' maxmiss90_common_snps.ped
sed -i 's/^SI_FIO[0-9]*/Fiordland/' maxmiss90_common_snps.ped
sed -i 's/^SI_WES[0-9]*/Westland/' maxmiss90_common_snps.ped
sed -i 's/^SI_NEL[0-9]*/Nelson/' maxmiss90_common_snps.ped
```

###### 13.12.18
There is an issue with those input files, because admixture throws an error (a very uninformative one) and exits. In the past I have prepared the plink file with stacks populations, so I tried that in boros:
```bash
module load stacks
populations -V maxmiss90_common_snps_HWE_LD.recode.vcf -O . -M population.txt --plink -t 8
mv maxmiss90_common_snps_HWE_LD.recode.p.plink.map maxmiss90_common_snps.map
mv maxmiss90_common_snps_HWE_LD.recode.p.plink.ped maxmiss90_common_snps.ped
```
But admixture does not seem to like these files either. From a quick read on some help pages, recoding the files in plink should help, but plink is likely going to pain me with my chromosome names being weird. This might work:
```bash
/usr/local/plink-1.9/plink --file maxmiss90_common_snps --recode 12 --out maxmiss90_common_snps --allow-extra-chr
```
That fixed it, admixture starts fine. Just out of curiosity I tried using the same recoding command on the vcftools outputs as well and it worked fine. Since the stacks command had removed some loci before outputting the plink files, for reasons not clear to me...I am using the vcftools outputs instead.
Now I want to fix the admixture script so that it runs through the Ks but also so that it runs 5 times with different random seeds. I "picked" 5 numbers randomly generated online, so that each K can be run with the same random seed and I can compare the different runs a bit more easily. I also need to fix the log names so that they don't overwrite.
I realised after a quick test run that I also need a way to get the outputs in different folders for each seed or they will overwrite. So I added some directory nesting system to split the outputs from different seeds.
Moving over to boros and starting.
```bash
scp ../ModPop_repo/admixture.sh boros:/data/denise/ModPop_analysis/pop_structure
screen -S ADMIX
bash admixture.sh
```

###### 14.12.18
Admixture ran fine, best K is always 1 in all random seed replicates. There is some variation between replicates, but it is quite minimal. I will plot the first 4 Ks to check that it is the same as last time. (Well, not K=1 of course...that would be silly). Choosing the random seed that had the lowest overall CV values, #11. Later on I will plot the CV values from the different runs.
```bash
scp boros:/data/denise/ModPop_analysis/pop_structure/admix_11/*.Q ./admix_11
cp ../../Kaka_GBS/GBS_scripts/admixture_plots.R ../ModPop_repo/
```
Everything looks exactly like the first time.
I will come back to this and fix the script so that it also produces a CV plot and so that the plots are made easily for each K.
For now, on to the next analysis.

##### DAPC/ADEGENET
###### 14.12.18
I copied to the repo the script I used before to run adegenet for my GBS data. Quickly using it to plot some DAPC analysis.
I need a file in the PLINK `.raw` format, so I am producing that from the previous plink file that is still sitting in boros, then moving it over here.
```bash
/usr/local/plink-1.9/plink --file maxmiss90_common_snps --recode A --out maxmiss90_common_snps --allow-extra-chr
scp boros:/data/denise/ModPop_analysis/pop_structure/maxmiss90_common_snps.raw .
scp boros:/data/denise/ModPop_analysis/pop_structure/maxmiss90_common_snps.map .
```
###### 17.12.18
All the tests I ran are in a separate report, called `GBS_adegenet.Rmd`. Nothing particularly new, but I managed to output a phylogenetic tree which looks a bit interesting, because the SI samples seem to "envelop" the NI samples. Not unexpectedly, the Kapiti and Zealandia samples are very nested within the tree (not one within the other though). The other cool thing is that there is a discriminant between the NI and the SI, I extracted the snps that contribute to it and I am looking at where they end up in the annotation.

##### TREEMIX
###### 17.12.18
First of all, I need to install this. I decided to use it on boros, so I will try installing it there. I have verified that the two prerequisites mentioned in the manual (boost and GSL) are present in boros. I downloaded the tarball from `https://bitbucket.org/nygcresearch/treemix/downloads/`. Next, I moved it to boros and installed it:
```bash
scp -r /Users/denisemartini/Downloads/Treemix/treemix-1.13.tar.gz boros:/home/denise/
tar -xvf treemix-1.13.tar.gz
cd treemix-1.13/
./configure --prefix /home/denise  # specify prefix, otherwise it tries to put the executables in /usr/local and I don't have permission
make
make install
```
Then, on to prepare input files. They provide a python script with the distribution that converts from a plink frequencies format to the treemix input files. But I need to get that plink first. I need to fix the population file in a format that plink likes, as such:
```bash
# it needs to look like:
FAMILYID  INDVID  CLUSTER
# so I am adding a column before my individual IDs in the pop file, and I will need that to be the population name, because it is that way in the .ped file
cat population.txt | awk '{print $2, $1, $2}' > plink_clusters.txt
sed -i 's/Kapiti/KapitiIsland/g' plink_clusters.txt
sed -i 's/LittleBarrier/LittleBarrierIsland/g' plink_clusters.txt
sed -i 's/Codfish/CodfishIsland/g' plink_clusters.txt
```
Then the plink conversion should work:
```bash
/usr/local/plink-1.9/plink --file maxmiss90_common_snps --freq --missing --within plink_clusters.txt
# output is in plink.frq.strat
less plink.frq.strat
CHR                 SNP     CLST   A1   A2      MAF    MAC  NCHROBS
  0      un_sc_20349_15 CodfishIsland    1    2        0      0       10
  0      un_sc_20349_15 Fiordland    1    2    0.125      2       16
  0      un_sc_20349_15 KapitiIsland    1    2        0      0       26
  0      un_sc_20349_15 LittleBarrierIsland    1    2        0      0       24
  0      un_sc_20349_15   Nelson    1    2        0      0       14
  0      un_sc_20349_15  Pureora    1    2   0.1111      2       18
  0      un_sc_20349_15 Westland    1    2        0      0       26
  0      un_sc_20349_15 Zealandia    1    2        0      0       16
```
Looking okay. Then, following treemix instructions:
```bash
gzip plink.frq.strat
/home/denise/bin/plink2treemix.py plink.frq.strat.gz treemix.frq.gz
```
Seems to have run fine. Let's try with the f3 and f4 statistics now.
```bash
/home/denise/bin/threepop -i treemix.frq.gz -k 500 > f3_statistics.txt
/home/denise/bin/fourpop -i treemix.frq.gz -k 500 > f4_statistics.txt
```
Done. I don't see anything negative in the f3 statistics, but there is a lot different from 0 in the f4 statistics.
I am taking away some extra bits in between, moving the file back to HPC and filtering the significant values.
```bash
grep -v "Estimating" f4_statistics.txt > f4_statistics_fixed.txt
grep -v "total_nsnp" f4_statistics_fixed.txt > f4_statistics_fixed2.txt
scp boros:/data/denise/ModPop_analysis/pop_structure/f4_statistics_fixed2.txt .
awk -F' ' '($4>=3.29)' f4_statistics_fixed2.txt > f4_significant_05.txt
awk -F' ' '($4<=-3.29)' f4_statistics_fixed2.txt >> f4_significant_05.txt
awk -F' ' '($4>=3.71)' f4_statistics_fixed2.txt > f4_significant_01.txt
awk -F' ' '($4<=-3.71)' f4_statistics_fixed2.txt >> f4_significant_01.txt
```
MAYBE, after seeing these stats, the admixture is not so generalised that I can't see it in the treemix tree.
I am going to try.
```bash
/home/denise/bin/treemix -i treemix.frq.gz -o kaka_nomig
```
Then to visualise the tree and residues:
```R
source("/home/denise/treemix-1.13/src/plotting_funcs.R")
plot_tree("kaka_nomig")
plot_resid("kaka_nomig", "poporder")
```
I made the poporder file simply by listing the pops in the population file [`awk '{print $2}' population.txt | sort -u > poporder`].
There are discrepancies obviously, specifically between LittleBarrierIsland and Zealandia.
So, let's try adding some migration.
```bash
/home/denise/bin/treemix -i treemix.frq.gz -root LittleBarrierIsland -m 6 -o kaka_withmig
```
Plotting:
```R
source("/home/denise/treemix-1.13/src/plotting_funcs.R")
plot_tree("kaka_withmig")
plot_resid("kaka_withmig", "poporder")
```
This looks confusing, and it changes a lot if I move the root to a different outgroup. I need a proper outgroup, to do some bootstrap tests and to try removing some missing data. But the program is really fast! =)

##### Genetic diversity stats
###### 14.3.19
I just thought that it would be nice to have a table or something with the average heterozygosity and minor allele frequency (and nucleotide diversity?) per population. Maybe rates of shared or fixed polymorphic sites as well?
So, I will prepare a file with the individuals from each pop and use VCFtools to extract stats for each pop, then average them and put them all in a table.
```bash
mkdir gen_diversity
cd gen_diversity
pops="Kapiti LittleBarrier Pureora Zealandia Codfish Nelson Fiordland Westland"
for p in $pops
do
awk -v var="$p" '$0~var{print $1}' ../population.txt > $p.txt
done
```
Then to get stats:
```bash
cp ../maxmiss90_common_snps_HWE_LD.recode.vcf .
vcftools=/usr/local/vcftools_0.1.15/bin/vcftools # since I am working in boros
$vcftools --vcf maxmiss90_common_snps_HWE_LD.recode.vcf \
--het --out all_het # this is individual, I don't need to use the split pops for this
for p in $pops
do
$vcftools --vcf maxmiss90_common_snps_HWE_LD.recode.vcf \
--keep ${p}.txt --freq2 --out ${p}
done
```
I think this is all I need to do some manipulation in R. Log in `GBS_gendiversity_stats.Rmd`.
Moving the files to desktop so I can work in R over there.
`scp -r boros:/data/denise/ModPop_analysis/pop_structure/gen_diversity ./pop_structure`

##### Tree using SNAPP/BEAST
###### 21.3.19
I need to do some filtering on my SNPs, because SNAPP is super computationally intensive. I will have to limit my number of samples to 3 from each population (24 total) and might still be too many. I will decide which individuals to keep based on their levels of missing data. I can obtain this information from VCFtools.
```bash
mkdir SNAPP_tree
cd SNAPP_tree
VCFtools --vcf ../vcf_filtering/maxmiss90_common_snps_HWE_LD.recode.vcf \
--missing-indv --out missingness_indv
awk '{print $1,$5}' missingness_indv.imiss | sed 's/\([A-Z]*_[A-Z]*\)\([0-9]*\) /\1\2 \1 /' | sort -k2,3 > sorted_indv_missingness.txt
```
Deleting all but the first three individuals from each pop in nano. One of the Nelson individuals still has 1/3 missing sites unfortunately, but that will probably not be a problem when I am done reducing the dataset. For now, let's put these individuals in a file, filter everything else out with VCFtools and keep only what still is a polymorphic site.
```bash
awk '{print $1}' sorted_indv_missingness.txt > selected_indv.txt
VCFtools --vcf ../vcf_filtering/maxmiss90_common_snps_HWE_LD.recode.vcf \
--keep selected_indv.txt --maf 0.001 --recode --out selected_indv_snps
After filtering, kept 80317 out of a possible 101539 Sites
```
That is still way too many, so let's restrict the missing sites as well.
```bash
VCFtools --vcf ../vcf_filtering/maxmiss90_common_snps_HWE_LD.recode.vcf \
--keep selected_indv.txt --maf 0.001 --max-missing 1 --recode --out selected_indv_snps
After filtering, kept 36425 out of a possible 101539 Sites
```
Still a lot! I can probably increase the maf a bit still: with 24 individuals and 48 alleles, without missing data at all, if I want at least 4 alleles to be the minor allele, that's a maf of 0.083.
```bash
VCFtools --vcf ../vcf_filtering/maxmiss90_common_snps_HWE_LD.recode.vcf \
--keep selected_indv.txt --maf 0.08 --max-missing 1 --recode --out selected_indv_snps
After filtering, kept 13551 out of a possible 101539 Sites
```
Getting there, I guess...I really only want to keep ~5000 SNPs. Short of randomly sampling them, the only other thing I can think of is the site quality. But this would not be even across the different programs I used. So, let's go for random, but I will need to do some roundabout trick.
```bash
grep -v '#' selected_indv_snps.recode.vcf > snps_to_sample.vcf
wc -l snps_to_sample.vcf
13551
grep '#' selected_indv_snps.recode.vcf > vcf_header.txt
cat snps_to_sample.vcf | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .373) print $0}' > sampled_snps.txt
# I had to play and rerun the command above a few time to get close to a number I was happy with,
# since it randomly extracts a proportion of lines (<= .373), not a fixed number
wc -l sampled_snps.txt
5004
# I just removed 4 extra random lines to get to 5000 SNPs, then put the file back together:
cat vcf_header.txt sampled_snps.txt > final_sampled_snps.vcf
```
Now I can use VCFtools to convert this to a .ped file, which I can then make into a fasta alignment.
```bash
VCFtools --vcf final_sampled_snps.vcf \
--plink --out final_sampled_snps
cut -f 1,7- final_sampled_snps.ped | tr -d '\t' | sed 's/\([A-Z]*_[A-Z]*\)\([0-9][0-9]\)/>\1\2 /' | tr ' ' '\n' > final_snp_alignment.fasta
cat final_snp_alignment.fasta | awk 'NR%2==0' | awk '{print length($1)}' #sanity check
10000
10000
10000
```
Nope, I had not considered one important factor: ped files actually have 2 columns per snp, the full genotype. What I need is a binary file, with SNPs coded as 0,1,2s depending on number of non-ref alleles. I think I can use the VCF 012 output and work from there to get my nexus file directly?
```bash
VCFtools --vcf final_sampled_snps.vcf \
--012 --out final_sampled_snps
cat final_sampled_snps.012 | cut -f 2- | tr -d '\t' > fixed_final_sampled_snps.012
paste final_sampled_snps.012.indv fixed_final_sampled_snps.012 > final_snp_matrix.012
```
Then I am adding in nano the bits at the beginning and at the end that make it a nexus file:
```
#NEXUS
BEGIN Data;
DIMENSIONS NTAX=24 NCHAR=5000;
FORMAT DATATYPE=integer INTERLEAVE=no missing=?;
Matrix
...<here>...
...<goes>...
...<mydata>...
;
END;
```
Saving as final_snp_matrix.nex
Now I am moving this in boros, because I can start BEAST over there.
```bash
cd ..
scp -r SNAPP_tree boros:/data/denise/ModPop_analysis
module load beast/2.5.0
beauti
beast -seed 777 -beagle -beagle_CPU -beagle_double -window -threads 12
```
###### 21.3.19
Right, so that is not going to work. Well, it is, but in about two months, which is not ideal right now. Let's restrict the dataset still.
First, I am going back to that selected samples file and leaving only the first two, then repeating the following selections.
```bash
awk '{print $1}' sorted_indv_missingness.txt > selected_indv.txt
VCFtools --vcf ../vcf_filtering/maxmiss90_common_snps_HWE_LD.recode.vcf \
--keep selected_indv.txt --maf 0.1 --max-missing 1 --recode --out selected_indv_snps
After filtering, kept 11099 out of a possible 101539 Sites
```
Now picking only 2500 SNPs.
```bash
grep -v '#' selected_indv_snps.recode.vcf > snps_to_sample.vcf
wc -l snps_to_sample.vcf
11099
grep '#' selected_indv_snps.recode.vcf > vcf_header.txt
cat snps_to_sample.vcf | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .225) print $0}' > sampled_snps.txt
# I had to play and rerun the command above a few time to get close to a number I was happy with,
# since it randomly extracts a proportion of lines (<= .373), not a fixed number
wc -l sampled_snps.txt
2505
# I just removed 5 extra random lines to get to 5000 SNPs, then put the file back together:
cat vcf_header.txt sampled_snps.txt > final_sampled_snps.vcf
```
Finally, getting the right output from VCFtools.
```bash
VCFtools --vcf final_sampled_snps.vcf \
--012 --out final_sampled_snps
cat final_sampled_snps.012 | cut -f 2- | tr -d '\t' > fixed_final_sampled_snps.012
paste final_sampled_snps.012.indv fixed_final_sampled_snps.012 > final_snp_matrix.012
```
Fixed the file in the nexus format as above, then moved to boros. I noticed that BEAST is also installed in NeSI, so I will try using it over there, but I don't think I have access to the graphical interface over there, so I will still prepare the .xml file in beauti in boros and then move it to NeSI, after figuring out how much time I might need as well.
```bash
module load beast/2.5.0
beauti
sed -i 's/totalcount="4"/totalcount="3"/' final_snp_matrix.xml #beauti for some reason messes up with my diploid data, so I need to add this fix before running beast
beast -seed 777 -beagle -beagle_CPU -beagle_double -window -threads 12
```
I am running it briefly, just to see how long it still takes to do a series of samples, so I can maybe estimate how to use it in NeSI. The good thing is that BEAST does save the states as it goes, so if it gets interrupted it doesn't have to start from the beginning I think.
While that runs I will see how I can make this work in NeSI.
```bash
module load BEAST/2.4.7
beast
```
Turns out, that when I typed beast I did get the GUI popping up. Which is fantastic, because it means I don't have to go crazy with figuring out how to install SNAPP without beauti. So, I did that. It would still be less than optimal to use the GUI in NeSI, so I will forget about that and remove the "-window" option from my beast command. I will put this in a NeSI script:
```
#!/bin/bash -e
#SBATCH --job-name=BEAST      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=24:00:00         # Walltime (HH:MM:SS)
#SBATCH --cpus-per-task=24      # number of cores per task
#SBATCH --mem-per-cpu=1500      # memory/cpu (in MB)
#SBATCH --ntasks=2              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis/SNAPP_tree   # directory where you run the job
#SBATCH --output=%x-%j.out      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

module load BEAST/2.4.7
beast -seed 777 -beagle -beagle_CPU -beagle_double -threads 48 final_snp_matrix.xml
```
And let's just see what happens when I run this.
```bash
sbatch nesi_beast.sh
Submitted batch job 2877945
```
Because in the meantime, on boros, beast ran 1500 samples in ~20 minutes. Which I believe would still take ~15 days.
Comparatively, in the last ~5 minutes on NeSI beast has already run 2700 samples. I am optimistic. I think it will still take about 2 days, but that's still a huge improvement. Hopefully the chain is long enough to reach convergence.

###### 26.3.19
Right, so, as expected one day was not enough to run the full chain, so I had resumed the run by changing the nesi command to:
```bash
beast -seed 777 -beagle -beagle_CPU -beagle_double -threads 48 -resume final_snp_matrix.xml
```
But turns out that when you do that it assumes you just want more chains added and it ignores the previous target. Or I am assuming that that's what happened, because it just kept going...I even restarted the job a third time, and it still kept going...so I stopped it when I saw it had reached about 2.5M chains. Also because I had in the meantime taken a look at the traces and noticed that it had already reached convergence and sufficient ESS values. I installed and launched tracer like so:
```bash
wget https://github.com/beast-dev/tracer/releases/download/v1.7.1/Tracer_v1.7.1.tgz
tar -xzvf Tracer_v1.7.1.tgz
chmod -R 777 Tracer_v1.7.1
rm *.tgz
/nesi/project/uoo02327/programs/Tracer_v1.7.1/bin/tracer
module load BEAST/2.4.7
densitree   # and made a tree
```
To really ensure convergence, I prepared in beauti another run exactly like the first, but with a different output file, and I am running it like:
```bash
beast -seed 345129842 -beagle -beagle_CPU -beagle_double -threads 48 final_snp_matrix_chain2.xml
```
This time I am also giving full 48hrs to the script, which should be enough to run the full chain the first time.

###### 28.3.19
I prepared all the other xml files I need to run with the taxa assigned by population and by island.
Running the first chain of the pop taxa:
```
beast -seed 8660864 -beagle -beagle_CPU -beagle_double -threads 48 final_snp_matrix_pop_chain1.xml
```
I am giving it 50 hrs, because it looks like it does not quite get to the end of the chain safely in 48hrs.
The commands for the rest of the runs will be:
```
beast -seed 1299718 -beagle -beagle_CPU -beagle_double -threads 48 final_snp_matrix_pop_chain2.xml
beast -seed 4722735 -beagle -beagle_CPU -beagle_double -threads 48 final_snp_matrix_isl_chain1.xml
beast -seed 2244459 -beagle -beagle_CPU -beagle_double -threads 48 final_snp_matrix_isl_chain2.xml
```
The seeds are randomly generated, of course.

##### Migration surfaces with EEMS
###### 25.3.19
EEMS is a pretty cool program that I first heard about at SMBE last year, and the name stands for "estimated effective migration surface" https://github.com/dipetkov/eems. What it does is take SNP information and geographical distribution of samples and plot on a map the migration surfaces for that population/species, allowing you to identify possible barriers to migrations and in general areas within there is more or less gene flow than expected. I am keen to try it out, and I would expect to find no substantial barriers to migration in my dataset. If the pattern is of Isolation-By-Distance, there shouldn't be any real barrier to migration, and basically no migration surfaces plottable because every area interacts with all others just as expected by sheer geography.
I have already checked that boros has the two dependencies (Eigen and Boost) installed and the versions should be fine. I need to download the code from github and compile it.
```bash
wget https://github.com/dipetkov/eems/archive/master.zip
unzip master.zip
cd eems-master/runeems_snps/src/
make linux
```
That failed, and it seems to be a problem with boost version. So, I saw suggestion online to install EEMS within a conda environment with the correct library versions.
```bash
module load conda
conda create -n eems
source activate eems
conda install boost=1.57
conda install eigen
cd eems-master/runeems_snps/src/
nano Makefile # then substitute these bits at the beginning of the file:
EIGEN_INC = /home/denise/.conda/envs/eems/include/eigen3
BOOST_LIB = /home/denise/.conda/envs/eems/lib
BOOST_INC = /home/denise/.conda/envs/eems/include
# save and
make linux
```
Then, to install the small program that comes with it and that should be used to covert plink files to the format they want:
```bash
cd ~
git clone https://github.com/mfranberg/libplinkio
cd libplinkio
git checkout 781e9ee37076
pip install plinkio
mkdir build
cd build
../configure --prefix=/home/denise/bin/plinkio
make && make check && make install
gcc -lplinkio -I/home/denise/bin/plinkio/include -L/home/denise/bin/plinkio/lib source.c
cd ../../
cd eems-master/bed2diffs/src-wout-openmp/
nano Makefile # substitute these lines in the file:
PLINKIO = /home/denise/bin/plinkio
CXXFLAGS = -I${PLINKIO}/include -O3 -Wall -Werror
# save and run
make linux
cd ../src/
nano Makefile # substitute these lines in the file:
PLINKIO = /home/denise/bin/plinkio
CXXFLAGS = -I${PLINKIO}/include -O3 -Wall -Werror
# save and run
make linux
```
Now, to prepare the input file, I first need to convert my dataset to plink binary file. I will use as input the finished plink .ped files that I produced before for admixture.
```bash
cd /data/denise/ModPop_analysis/pop_structure
mkdir EEMS
cd EEMS
/usr/local/plink-1.9/plink --file ../maxmiss90_common_snps --make-bed --out eems_snp_dataset
# now running the eems conversion
/home/denise/eems-master/bed2diffs/src/bed2diffs_v1 --bfile eems_snp_dataset
```
Seems to have worked fine.
Next, I need a file with the coordinates for each sample, one sample per line.
I prepared a similar file for the environmental correlations a while ago. I need to sort the samples by name, because that was the order of individuals in the .bed file that produced the snp dataset. I also need to have only the longitude and latitude of each sample.
```bash
cp ../../env_correlations_maf05/samples_locations.txt .
tail -n +2 samples_locations.txt | sort | awk -F'\t' '{print $4,$3}' | sed -e "s/\r/\t/g" > eems_snp_dataset.coord
```
Now I only need the outer coordinates file, which I need to put together. It needs to have the coordinates of the whole habitat circle, listed counterclockwise and in a closed circle, with the last coordinate being the same as the first.
I used this online tool https://www.keene.edu/campus/maps/tool/ to get the coordinates surrounding New Zealand (N.B. make sure to position the map correctly, from the US starting point and rotating towards the right, not the left of the map to get to NZ). I decided to include the whole of New Zealand, since that was the kaka distribution until not too long ago. I put these coordinates in a file and fixed them (to remove commas):
```bash
nano eems_snp_dataset.outer
sed -i 's/,//g' eems_snp_dataset.outer
```
Now everything should be in place.
Then I need to prepare three parameter files to run the mcmc chain independently three times. I can copy the files from the eems distribution:
```bash
cp /home/denise/eems-master/runeems_snps/src/params-chain*.ini .
nano params-chain1.ini # and filling with the following parameters:
datapath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset
mcmcpath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset-chain1
nIndiv = 92
nSites = 101539
nDemes = 200
diploid = true
numMCMCIter = 2000000
numBurnIter = 1000000
numThinIter = 9999
```
I did the same for the other two parameter files, just changing the number of the chain in the output file name ("mcmcpath").
Now, to run it:
```bash
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain1.ini --seed 123
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain2.ini --seed 456
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain3.ini --seed 789
```
I started all three iterations in separate tmux windows. It seems to be running pretty fast but it has quite a lot of iterations to go through, so I don't expect it to finish before tomorrow. But no errors at least.

###### 26.3.19
It actually finished running in an hour or so, way faster than I expected. So, I will need to rerun with slightly different parameters, because the first time they were accepted too often (see manual, section 2.1), but first I want to see if the length of the chain was enough to reach convergence, so I need the plotting functions that come from their R package. I need to install a couple of dependencies as well.
```R
install.packages("RcppEigen")
install.packages("rgeos", configure.args = c(
                    "--with-geos-config=/home/denise/bin/bin/geos-config"))
setwd("~/eems-master/plotting")
if (dir.exists("rEEMSplots")) {
  install.packages("rEEMSplots", repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}
```
Then to produce the plot results:
```R
setwd("/data/denise/ModPop_analysis/pop_structure/EEMS")
library("rEEMSplots")
eems.plots(mcmcpath = c("eems_snp_dataset-chain1", "eems_snp_dataset-chain2", "eems_snp_dataset-chain3"),
  plotpath = "eems_results-default",
  longlat = TRUE)
```
Looking cool and interesting so far, but also like it has not really converged yet, so I will fix parameters and rerun the whole thing. Since it is so fast I can afford to lengthen the chains quite a bit.
```
nano params-chain1.ini # and filling with the following parameters:
datapath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset
mcmcpath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset-chain1
nIndiv = 92
nSites = 101539
nDemes = 400
diploid = true
numMCMCIter = 5000000
numBurnIter = 1000000
numThinIter = 9999
qEffctProposalS2 = 0.002
mEffctProposalS2 = 0.15
mSeedsProposalS2 = 0.02
```
I also increased the number of demes, and that might increase the computational costs quite a bit. I fixed all three parameter files and let's run:
```bash
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain1.ini --seed 123
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain2.ini --seed 456
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain3.ini --seed 789
```
###### 27.3.19
Still not quite right, the acceptance rate is still a bit too high, so:
```
nano params-chain1.ini # and filling with the following parameters:
datapath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset
mcmcpath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset-chain1
nIndiv = 92
nSites = 101539
nDemes = 400
diploid = true
numMCMCIter = 5000000
numBurnIter = 1000000
numThinIter = 9999
qEffctProposalS2 = 0.004
mEffctProposalS2 = 0.25
mSeedsProposalS2 = 0.025
```
And restarting it as usual. If it runs like last night it should take ~6hrs.

###### 28.3.19
One last time:
```
nano params-chain1.ini # and filling with the following parameters:
datapath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset
mcmcpath = /data/denise/ModPop_analysis/pop_structure/EEMS/eems_snp_dataset-chain1
nIndiv = 92
nSites = 101539
nDemes = 400
diploid = true
numMCMCIter = 5000000
numBurnIter = 1000000
numThinIter = 9999
qEffctProposalS2 = 0.006
mEffctProposalS2 = 0.6
mSeedsProposalS2 = 0.035
```
And restarting it again.

###### 29.3.19
Parameters are not quite perfect yet, but much closer to where they should be. Let's take another look at those plots.
```R
library("rEEMSplots")
library("rgdal")
library("rworldmap")
library("rworldxtra")
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
eems.plots(mcmcpath = c("eems_snp_dataset-chain1", "eems_snp_dataset-chain2", "eems_snp_dataset-chain3"),
  plotpath = "eems_results-default",
  longlat = TRUE, add.grid = TRUE, col.grid = "gray90", lwd.grid = 0.8,
  add.demes = TRUE, col.demes = "black", pch.demes = 5, min.cex.demes = 0.5, max.cex.demes = 1.5,
  projection.in = projection_none, projection.out = projection_mercator, add.map = TRUE, col.map = "black", lwd.map = 1)
```
###### 1.4.19
After a few other tweaks over the weekend, this is I think the winning combination of parameters:
```
qEffctProposalS2 = 0.008
mEffctProposalS2 = 0.95
mSeedsProposalS2 = 0.04
```
###### 5.4.19
Those are winning parameters I think, but looking at the plots again, maybe convergence has not quite been reached yet. I will run it again later with a longer chain. But first, Michael and I wanted to see what would happen to these migration surfaces if Kapiti and Zealandia were not there, since they seem to be a hole of non-migration. So, I need to take them out from the input files.
```bash
mkdir noKZ
cp eems_snp_dataset.outer noKZ/ #this stays the same
cp eems_snp_dataset.coord noKZ/ #from this I remove Kapiti and Zealandia samples in nano
cp ../maxmiss90_common_snps.ped noKZ/ #also removed Kapiti and Zealandia in nano
cp ../maxmiss90_common_snps.map noKZ/ #this stays the same
/usr/local/plink-1.9/plink --file maxmiss90_common_snps --make-bed --out eems_snp_dataset
# now running the eems conversion
/home/denise/eems-master/bed2diffs/src/bed2diffs_v1 --bfile eems_snp_dataset
cp ../params-chain*.ini .
```
Fixing just the first lines of the params files:
```
datapath = /data/denise/ModPop_analysis/pop_structure/EEMS/noKZ/eems_snp_dataset
mcmcpath = /data/denise/ModPop_analysis/pop_structure/EEMS/noKZ/eems_snp_dataset-chain1
nIndiv = 67
```
I also shortened the chain a bit, I don't want it to run forever. I left all other params to the latest optimization ones.
```
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain1.ini --seed 123
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain2.ini --seed 456
/home/denise/eems-master/runeems_snps/src/runeems_snps --params params-chain3.ini --seed 789
```
The results are pretty neat. I mean, they haven't quite converged yet, which is expected, but the real nice thing is that I was right, if you exclude Kapiti and Zelandia the Cook strait does not look like a barrier to migration anymore. There is still the Nelson population looking grim in there, but I did have lower diversity there as well in my stats. So, I think I will try removing that as well, and my guess is that then I will see pretty much an IBD pattern. I need to remember that the Nelson population has also been probably isolated for a while, since the current kaka distribution pretty much excludes the central SI. So, the migration surfaces are pretty much reflecting the current distribution and its fragmentations. I will try that out, to confirm my interpretation. Running everything just like before, but after removing the Nelson samples from the input files and increasing the chains (it can run overnight).
```
nIndiv = 59
numMCMCIter = 8000000
```
Those were the only two things I changed in the params files.

###### 8.4.19
I am rerunning the normal dataset (with all individuals) for a bit longer, increasing the chain length to 8000000. Left all other parames as before.

##### Modeling in dadi
###### 26.3.19
I have been studying this for a while now and I have been wanting to try out this program ever since I found out about it in 2016. It is not an easy one, but might be one of the most interesting to use with this dataset. So, there is some preparation to do. dadi is basically a python language, so to run it you write a script of the model you want to test. Luckily, I found a nice resource online where all the scripts for the possible modelsI would like to test have already been implemented, together with a few other wrapping options: https://github.com/dportik/dadi_pipeline
I also found here and there some perl scripts that I should be able to use to go from my vcf files to the input required by dadi. The only thing I want to fix before I dig into it, is that I would like to have the option of using the information from the outgroup, the kea in my case, in the spectrum. Ergo, I need to call these SNPs in kea. I think I can do that from the whole-genome alignment of kea that I did for the species comparison analysis. I have everything in NeSI.
```bash
cd pop_structure
mkdir dadi
cp maxmiss90_common_snps_HWE_LD.recode.vcf dadi
cd dadi
# just to check what the chromosome names where in the original alignment
module load SAMtools
samtools idxstats ../../../Kea-Kaka_take2/realignment/kea_sorted_rmdup.bam > kea_idxstats.txt
# to call the positions from the vcf file (and fix chromosome names):
module load BCFtools
grep -v '#' maxmiss90_common_snps_HWE_LD.recode.vcf | awk '{print $1,$2}' | tr ' ' '\t' | sed 's/ch/chr/' > targets.txt
# then to call the genotypes at those positions:
samtools mpileup --positions targets.txt --BCF \
--uncompressed --output-tags DP,AD,INFO/AD \
--fasta-ref ../../../Kea-Kaka_take2/realignment/N_meridionalis_pseudochr.fa \
../../../Kea-Kaka_take2/realignment/kea_sorted_rmdup.bam | bcftools call \
--format-fields GQ,GP \
--multiallelic-caller \
--output-type v - --output kea_GBS_snps.vcf
```
Now, then, to add my kea vcf to my kaka vcf. First I need to bgzip and tabix both files, then merge them.
```bash
# bit of clean up on the kea file
sed -i 's/ps_chr/ps_ch/' kea_GBS_snps.vcf
sed -i 's/SM:kea/kea/' kea_GBS_snps.vcf
# getting rid of all unnecessary information from the vcfs:
bcftools annotate -x INFO,^FORMAT/GT kea_GBS_snps.vcf > clean_kea_GBS_snps.vcf
bcftools annotate -x INFO,^FORMAT/GT maxmiss90_common_snps_HWE_LD.recode.vcf > clean_kaka_GBS_snps.vcf
# bgzip and tabix
bgzip clean_kea_GBS_snps.vcf
tabix -p vcf clean_kea_GBS_snps.vcf.gz
bgzip clean_kaka_GBS_snps.vcf
tabix -p vcf clean_kaka_GBS_snps.vcf.gz
# forgot that I don't need the indels in the kea, so getting rid of them now:
module load VCFtools
vcftools --remove-indels --vcf clean_kea_GBS_snps.vcf --remove-filtered-all --recode --out noind_clean_kea_GBS_snps.vcf
# I worked on this last file in R a bit, because there were some problems with the reference sequence in places, so I opened it up in R, checked the places with problems against the kaka snp file, fixed the ref and genotypes and then put everything back together
grep -v '##' clean_kaka_GBS_snps.vcf > noheader_kaka_GBS_snps.vcf
grep -v '##' noind_clean_kea_GBS_snps.vcf.recode.vcf > noheader_kea_GBS_snps.vcf #these I imported in R
grep '#' noind_clean_kea_GBS_snps.vcf.recode.vcf > kea_header
cat kea_header fixed_kea_GBS_snps > fixed_kea_GBS_snps.vcf #the fixed file I exported from R
###
bgzip fixed_kea_GBS_snps.vcf
tabix -p vcf fixed_kea_GBS_snps.vcf.gz
bcftools merge -m snps clean_kaka_GBS_snps.vcf.gz fixed_kea_GBS_snps.vcf.gz > snps_with_outgroup.vcf
```
That finally went fine. The R code I used to fix the kea refs problem is in `fixing_kea_refs.R`.

###### 28.3.19
Now, I need to do a couple of quick filterings on this set. Specifically, I need to remove any sites that are not not biallelic anymore and I might as well remove sites for which kea has missing data, since both of these will be excluded by dadi and I think they would mess up the file conversion step.
```bash
module load VCFtools
vcftools --vcf snps_with_outgroup.vcf --min-alleles 2 --max-alleles 2 \
--remove-filtered-all --recode --out biallelic_snps_with_outgroup
After filtering, kept 101317 out of a possible 101553 Sites
grep -v '#' biallelic_snps_with_outgroup.recode.vcf | awk -F'\t' 'BEGIN{OFS="\t";} {if ($102 == "./.") print $1,$2}' > missing_in_kea.txt
wc -l missing_in_kea.txt
2130
vcftools --vcf biallelic_snps_with_outgroup.recode.vcf \
--exclude-positions missing_in_kea.txt \
--remove-filtered-all --recode --out final_snps_for_dadi
After filtering, kept 99187 out of a possible 101317 Sites
```
Finally, to see whether the scripts I got from the dadi-user group work. The one to convert vcf files to dadi input comes from this thread https://groups.google.com/forum/#!searchin/dadi-user/vcf|sort:date/dadi-user/p1WvTKRI9_0/tqQQUJX7AwAJ, and this link: https://github.com/wk8910/bio_tools/blob/master/01.dadi_fsc/00.convertWithFSC/convert_vcf_to_dadi_input.pl. It should work with just the vcf file and a list of samples and populations.
```bash
awk '{print $1,$3}' ../population.txt | tr ' ' '\t' > dadi_poplist.txt
```
I added to that one line: `kea  out` to indicate that the kea sample is the outgroup. Now let's see if the script works.
```
perl convert_vcf_to_dadi_input.pl final_snps_for_dadi.recode.vcf dadi_poplist.txt
```
That worked like a charm, I am very grateful to whoever wrote it! I only fixed a typo in the header (Position) and changed the flanking positions in the sequence to be "-", as reported in dadi manual, rather than Ns.
The output file looks like:
```
Ref     OUT     Allele1 North   South   Allele2 North   South   Gene    Position        
-C-     -C-     C       79      82      T       11      4       ps_ch_1 9925    
-C-     -C-     C       93      62      G       5       6       ps_ch_1 10540   
-A-     -A-     A       15      5       G       3       5       ps_ch_1 23215   
-T-     -T-     T       96      74      C       2       12      ps_ch_1 29295   
-C-     -C-     C       88      74      T       10      12      ps_ch_1 35815   
-T-     -T-     T       10      9       A       0       1       ps_ch_1 39436   
-G-     -G-     G       90      71      T       6       5       ps_ch_1 47471   
-T-     -T-     T       96      78      C       2       4       ps_ch_1 60927   
```
Time to get dadi setup. I will prepare it as a conda environment, so that I don't mess up my other python environments. I need to check what version of python is needed. --> Yup, as expected, Python >2.5 is required, but not Python3. Very annoying.
```bash
conda create -n dadi python=2.7
source activate dadi
conda install numpy scipy matplotlib
conda install ipython
conda install -c bioconda dadi
conda install ipykernel # these last two I only need to run atom with this env
conda install -c conda-forge python-language-server
```
Hopefully the conda install is solid, because it was an incredibly easy install like this.
I will try running some tests and see what happens. I am playing around with it in the script `figuring_out_dadi.py`. I tried a couple of the examples that come with the dadi distribution, then I will try loading up my data and figuring out how to project down my snps.
Now that I have done that and read up some more, I want to try fitting my two populations separately at first, to see what kind of growth or not growth I can expect from them. I think it might be safe to try fitting these "simple" models here on my desktop.
So, I downloaded the repo with the optimization scripts from https://github.com/dportik/dadi_pipeline and I am going to use his script with the models form 1D spectrum that come with the dadi distribution. To do these things I will need to modify the script accordingly.
```bash
git clone https://github.com/dportik/dadi_pipeline.git # in my ~ folder
cd /Volumes/Denise/ModPop_analysis/pop_structure/dadi/
cp ~/dadi_pipeline/Optimize_Functions.py .
cp ~/dadi_pipeline/dadi_Run_Optimizations.py .
cp ~/gutenkunstlab-dadi-f2f4b565089a/dadi/Demographics1D.py .
```
The second script is the one that needs to be modified. Specifically, I am removing a bunch of writing (instructions and examples), setting up my run, specifying all the models I want to try fitting and optimizing them one at a time. The models I want to try are all in Demographics1D.py, which I am importing in the script, it should be callable from dadi itself, but I copied it here just in case.
```
# in the dadi environment
python dadi_Run_Optimizations_North.py | tee logNorth.txt
python dadi_Run_Optimizations_South.py | tee logSouth.txt
```
Which such few parameters it really seems to be going pretty fast, but I can see how I would like to run it on the cluster with more complex models.

###### 1.4.19
While that really was quick (about 1hr to run all models) I didn't really have an opportunity to look at the results properly. So, I checked the likelihoods and all that by running the summarizing script from the dadi pipeline:
```bash
mkdir North
mv North.* North
python Summarize_Outputs.py ./North
mkdir South
mv South.* South
python Summarize_Outputs.py ./South
```
So I get that for the NI the best model is the bottlegrowth, but with growth followed by bottleneck, while for the SI I get a three_epoch model, with prolonged bottleneck and recovery. I am not sure of how comparable the parameters for the two populations are...like the time or pop size parameters. But in both cases it seems like convergence has not quite been reached, so I will need to run this again, with more iterations or rounds of optimization. And the bottlegrowth model estimates for the popsize reached the upper bound parameters in both pops, so I will need to change those as well. But first I am taking a quick look at the simulations generated from these models, with the optimized parameters, to compare them with the real data and see how well (or not) they fit it.
I checked to fit of the models and I have a problem, probably due to ancestral misidentification, which seems to be fairly common with GBS/RAD datasets. So, following the suggestions here https://bit.ly/2Ul83f3 and here https://bit.ly/2FNyodA I am integrating a correction in my models. I will add a line in the dadi_Run_Optimizations scripts, that corrects the model function before running the optimizations and adds this parameter. I will need to remember that there is this extra parameter when I specify the optimization parameters.
```py
snps = "../final_snps_for_dadi.recode.vcf.data"
prefix = "North_misid"
pts = 60,70,80]
reps = [20,30,50]
maxiters = [10,15,25]
folds = [3,2,1]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.two_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "two_epoch", func_anc, 3, 3, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T, p_misid")
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.growth)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "growth", func_anc, 3, 3, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T, p_misid")
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.three_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch", func_anc, 3, 5, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF, p_misid")
upper = [50,30,30,30]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.bottlegrowth)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", func_anc, 3, 4, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, T, p_misid", in_upper=upper)
```
I am also increasing the number of repetitions and iterations for each round of optimizations. And increasing upper bounds for the bottlegrowth model.
There is just one problem: the conda installation of dadi does not include this function (it came with a later version of the program I guess). So let's fix that:
```bash
conda remove dadi
cd ~
# I went crazy because the package is in the middle of switching to Python3...and I had to find the specific commit where the developers started the switch -.-
wget https://bitbucket.org/gutenkunstlab/dadi/get/f2f4b565089af660573339d88efe9548d2a596fa.zip
unzip f2f4b565089af660573339d88efe9548d2a596fa.zip
mv gutenkunstlab-dadi-f2f4b565089a/ dadi
rm f2f4b565089af660573339d88efe9548d2a596fa.zip
cd dadi
python setup.py build_ext --inplace
# then I need to add this to the PYTHONPATH
export PYTHONPATH=$PYTHONPATH:/home/denise/dadi
```
Now, it should work:
```
cd /data/denise/ModPop_analysis/pop_structure/dadi
cd North
python -u ../dadi_Run_Optimizations_North.py | tee ../logNorth_misid.txt
cd ../South
python -u ../dadi_Run_Optimizations_South.py | tee ../logSouth_misid.txt
```
Then collecting results:
```
python Summarize_Outputs.py ./South
python Summarize_Outputs.py ./North
```
Still not quite there, the fit is close but still deviating from the real data at the derived and ancestral low frequency alleles. It might still be that the ancestral misidentification was not corrected enough (even though the parameter was optimized at ~4.3% in both populations), or that the nuB parameter (pop size after the bottleneck, or in this case after the rapid expansion) is still hitting the upper boundaries.
Even after increasing the boundaries again, I don't see a real improvement. Time to fold the SFS. I will fix the Optimization script again, to use the folded SFS and forget about the ancestral misidentification.
```py
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = False)
prefix = "North_folded"
pts = [60,70,80]
reps = [20,30,50]
maxiters = [10,15,25]
folds = [3,2,1]
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "two_epoch", Demographics1D.two_epoch, 3, 2, fs_folded=True, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T")
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "growth", Demographics1D.growth, 3, 2, fs_folded=True, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T")
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch", Demographics1D.three_epoch, 3, 4, fs_folded=True, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF")
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", Demographics1D.bottlegrowth, 3, 3, fs_folded=True, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, T")
```
###### 2.4.19
So, after folding the spectrum last night and rerunning, I feel like we are getting really close, the fit is almost perfect for the three epoch model in both populations. But it still looks like the optimization isn't quite there yet. And the results are still similar to the ones with the unfolded spectrum: very big growth followed by decline/bottleneck. I really wonder if I can't just get a good fit for the folded spectrum if I run the optimizations like, forever. Maybe starting with initial parameters close to the ones that were found in the best rounds last night?
```py
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections = proj, polarized = True)
prefix = "North_misid_long" # or "South_misid_long"
pts = [80,90,100]
reps = [20,30,40,50]
maxiters = [10,15,20,25]
folds = [3,2,2,1]

params = [2.8,1.2,0.04]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.two_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "two_epoch", func_anc, 4, 3, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T, p_misid", in_params=params)
params = [8.0,12.0,0.04]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.growth)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "growth", func_anc, 4, 3, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nu, T, p_misid", in_params=params)
params = [12.0,1.2,0.7,0.1,0.04]
upper = [100,30,30,30,30]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.three_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch", func_anc, 4, 5, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF, p_misid", in_params=params)
params = [75.0,1.3,0.8,0.04]
upper = [120,30,30,30]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.bottlegrowth)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "bottlegrowth", func_anc, 4, 4, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, T, p_misid", in_upper=upper, in_params=params)

# starting params for SI
params = [3.0,0.8,0.04] #two_epoch
params = [4.0,1.9,0.04] #growth
params = [9.3,0.7,0.5,0.1,0.04] #three_epoch
params = [35.0,2.0,0.8,0.04] #bottlegrowth
```

I think that took us quite close to the optimal fit. This is probably the best model I can get for my 1D populations, and it is relatively simple. But before diving into 2D models, I want to check a couple of other possibilities, adding some complexity to the three_epoch model in particular to see if it adds anything to it or not. Specifically, I am adding a couple more size changes, to verify if we are picking up a single recent bottleneck or two, and exponential rather than instantaneous size changes. These are several more models and I will test them all together against the three_epoch that we have already "established".
```py
import Exp1DModels
prefix = "North_ext1D"
params = [27.0,0.1,0.6,0.01,0.04]
lower = [0.01,0.01,0.01,0.001,0.01]
upper = [100,30,30,30,30]
func_anc = dadi.Numerics.make_anc_state_misid_func(Demographics1D.three_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch", func_anc, 4, 5, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF, p_misid", in_params=params, in_upper=upper, in_lower=lower)
func_anc = dadi.Numerics.make_anc_state_misid_func(Exp1DModels.four_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "four_epoch", func_anc, 4, 7, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuA, nuB, nuC, TA, TB, TC, p_misid")
func_anc = dadi.Numerics.make_anc_state_misid_func(Exp1DModels.five_epoch)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "five_epoch", func_anc, 4, 9, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuA, nuB, nuC, nuD, TA, TB, TC, TD, p_misid")
func_anc = dadi.Numerics.make_anc_state_misid_func(Exp1DModels.three_epoch_exp)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_exp", func_anc, 4, 5, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF, p_misid")
func_anc = dadi.Numerics.make_anc_state_misid_func(Exp1DModels.three_epoch_firstexp)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_firstexp", func_anc, 4, 5, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF, p_misid")
func_anc = dadi.Numerics.make_anc_state_misid_func(Exp1DModels.three_epoch_scndexp)
Optimize_Functions.Optimize_Routine(fs, pts, prefix, "three_epoch_scndexp", func_anc, 4, 5, fs_folded=False, reps = reps, maxiters = maxiters, folds = folds, param_labels = "nuB, nuF, TB, TF, p_misid")
```
I put the above in new optimization scripts, to avoid overwriting the last ones, which worked fairly well. Because these models are a bit more complex (especially the four and five epochs) I don't necessarily expect them to converge right away, but this might be a start.
```
scp ../../ModPop_repo/Exp1DModels.py boros:/data/denise/ModPop_analysis/pop_structure/dadi
python -u ../dadi_Run_Ext1D_NI.py | tee ../logNorth_misid.txt
python -u ../dadi_Run_Ext1D_SI.py | tee ../logSouth_misid.txt
```
From a quick look at results I can say that it's probably not worth going to the effort of trying to optimize the four and five_epoch models, too many parameters and there isn't much of an increase in likelihood I could gain at this stage. The three epoch models with exponential size changes instead might be a good option, they get close even if they are not optimized properly. So, I might consider including some of that in the models I test with the 2D spectrum. I will put a few of the most basic models in a script, then run veeeeery long optimizations with them.

I set up the models in basic_2DModels.py and I split the optimization routines in two scripts, to run them at the same time. I set them up to be fairly long and intensive and even increased the grid size. Some of these models are very parameter rich, so I expect this will take a while, especially on the 2D spectrum. Hopefully it will run fine.
```
scp ../../ModPop_repo/basic_2DModels.py boros:/data/denise/ModPop_analysis/pop_structure/dadi
scp ../../ModPop_repo/dadi_Run_2D_Set*.py boros:/data/denise/ModPop_analysis/pop_structure/dadi
mkdir 2Dfs
cd 2Dfs
python -u ../dadi_Run_2D_Set1.py | tee ../log2D_Set1.txt
python -u ../dadi_Run_2D_Set2.py | tee ../log2D_Set2.txt
```
###### 3.4.19
Since this setup looked like it was going to take forever and ever, I reduced quite a few parameters (smaller grid size, less rounds, less reps per round and narrower parameter boundaries), to try and get a quick and dirty look at what the models look like they are fitting before starting the real long optimizations. Hopefully this will speed up things in the long term. Also, I noticed that the log files were not picking up the errors/warnings from the run and I would really like to have those included, so fixing the logging method now. I hope.
```
scp ../../ModPop_repo/dadi_Run_2D_Set*.py boros:/data/denise/ModPop_analysis/pop_structure/dadi
python -u ../dadi_Run_2D_Set1.py 2>&1 | tee ../log2D_Set1.txt
python -u ../dadi_Run_2D_Set2.py 2>&1 | tee ../log2D_Set2.txt
```
###### 4.4.19
That went much better, with a lot less errors than the first time as well, even though clearly we did not reach convergence. I will run this another couple of times, taking care not to overwrite the results, so that I can see if the optimizations are at least all heading in the same direction. So the only things I am changing are the prefix for the output files in the optimization routines.
```
prefix = "2Dfs_2"
scp ../../ModPop_repo/dadi_Run_2D_Set*.py boros:/data/denise/ModPop_analysis/pop_structure/dadi
python -u ../dadi_Run_2D_Set1.py 2>&1 | tee ../log2D_2_Set1.txt
python -u ../dadi_Run_2D_Set2.py 2>&1 | tee ../log2D_2_Set2.txt
prefix = "2Dfs_3"
scp ../../ModPop_repo/dadi_Run_2D_Set*.py boros:/data/denise/ModPop_analysis/pop_structure/dadi
python -u ../dadi_Run_2D_Set1.py 2>&1 | tee ../log2D_3_Set1.txt
python -u ../dadi_Run_2D_Set2.py 2>&1 | tee ../log2D_3_Set2.txt
```
As for the actual serious optimizations, I think I might like to run that in NeSI instead, because it might take much longer than this. So, I will install the environment I need over there as well...
```bash
cd /nesi/project/uoo02327/programs
module load Miniconda3/4.4.10
conda create -p /nesi/project/uoo02327/programs/miniconda_envs/dadi python=2.7
source activate /nesi/project/uoo02327/programs/miniconda_envs/dadi
conda install numpy scipy matplotlib ipython
wget https://bitbucket.org/gutenkunstlab/dadi/get/f2f4b565089af660573339d88efe9548d2a596fa.zip
unzip f2f4b565089af660573339d88efe9548d2a596fa.zip
mv gutenkunstlab-dadi-f2f4b565089a/ dadi
rm f2f4b565089af660573339d88efe9548d2a596fa.zip
cd dadi
python setup.py build_ext --inplace
# then I need to add this to the PYTHONPATH when I am running scripts
export PYTHONPATH=$PYTHONPATH:/nesi/project/uoo02327/programs/dadi
```
Then I should be able to run my dadi scripts...I might test out the fast ones for a second, because it might be convenient to split the models in more separate jobs. I just don't know if NeSI will be faster if not by splitting the job, since dadi does not actually run in parallel, so I don't know that giving it more memory in NeSI will help.
```
scp final_snps_for_dadi.recode.vcf.data mahuika://nesi/nobackup/uoo02327/denise/ModPop_analysis/pop_structure/dadi
scp Optimize_Functions.py mahuika://nesi/nobackup/uoo02327/denise/ModPop_analysis/pop_structure/dadi
scp ../../ModPop_repo/dadi_Run_2D_Set*.py mahuika://nesi/nobackup/uoo02327/denise/ModPop_analysis/pop_structure/dadi
scp ../../ModPop_repo/basic_2DModels.py mahuika://nesi/nobackup/uoo02327/denise/ModPop_analysis/pop_structure/dadi
scp Summarize_Outputs.py mahuika://nesi/nobackup/uoo02327/denise/ModPop_analysis/pop_structure/dadi
```
After moving there the essentials, I can just briefly put the commands to run the optimizations in a nesi script.
```
#!/bin/bash -e
#SBATCH --job-name=dadi_2D_Set2      # job name (shows up in the queue)
#SBATCH --account=uoo02327     # Project Account
#SBATCH --time=02:00:00         # Walltime (HH:MM:SS)
#SBATCH --mem=12GB      # memory
#SBATCH --ntasks=1              # number of tasks (e.g. MPI)
#SBATCH --partition=large       # specify a partition
#SBATCH --hint=nomultithread    # don't use hyperthreading
#SBATCH --chdir=/nesi/nobackup/uoo02327/denise/ModPop_analysis/pop_structure/dadi/2Dfs   # directory where you run the job
#SBATCH --output=%x-%j.log      # %x and %j are replaced by job name and ID
#SBATCH --mail-type=ALL         # Optional: Send email notifications
#SBATCH --mail-user=marde569@student.otago.ac.nz     # Use with --mail-type option

# loading conda and dadi
module load Miniconda3/4.4.10
source activate /nesi/project/uoo02327/programs/miniconda_envs/dadi
export PYTHONPATH=$PYTHONPATH:/nesi/project/uoo02327/programs/dadi

# running the optimization script
python -u ../dadi_Run_2D_Set2.py
```
I have tried asking directly for memory rather than cpu, not sure if that changes anything. 12GB of RAM should be about 8 CPUs in NeSI? I will leave this running for a little bit, then stop it, increase memory requirements and see if it gets any further...
I tested with 12GB and 24GB of mem requirements and it does not seem to make any difference. Seems to go smoother than in boros anyway, and to use far less memory than those requirements (probably ~250MB?). There is more variability in terms of time per model in boros, which makes me think that there might be some struggle there which might not be happening over here. The models simply don't seem to get stuck the way they do there. So, my final idea is that I will split my optimization runs a bit further, setting up the models in 4 jobs that should take about the same time: one with the easier SNM and two epoch models, one with the simple three epochs, one with the beforeafters and one with the three_epoch_afters. I will only ask for one CPU each and not too many hours for now (6?). And I will fix those optimization reps to be longer and with a larger grid size, as soon as I have an idea of what the best initial parameters might be.
```
scp ../../ModPop_repo/dadi_Run_2D_Set*-*.py mahuika://nesi/nobackup/uoo02327/denise/ModPop_analysis/pop_structure/dadi
sbatch nesi_dadi2D_Set1-1.sh
sbatch nesi_dadi2D_Set1-2.sh
sbatch nesi_dadi2D_Set2-1.sh
sbatch nesi_dadi2D_Set2-2.sh
```
###### 5.4.19
So, with slightly optimized parameters I seem to be getting less "model is masked" errors, but I think I am still getting lots (more?) "extrapolation failed" errors. Even though I did increase the grid size before running these. Also, the runs did not finish because there was the thing again where it got stuck on some repetitions. I will try switching to the linear extrapolation function, rather than the logarithmic one, by changing this in the Optimize_Functions.py script.
```py
# at line 233 change from:
func_exec = dadi.Numerics.make_extrap_log_func(func)
#to:
func_exec = dadi.Numerics.make_extrap_func(func)
```
I will also increase the time requirements a bit, just in case...I am asking for 10 hrs now. That problem of some models getting stuck seems really unpredictable though.
###### 8.4.19
After a few other trials over the last couple of days, I had to resort to put some of the models that were getting "stuck" more frequently in their own run, and asked for full 24hrs for each of them. It was a good idea, since one of the more complex models still managed to run for 21.5hrs on its own. But, even with that, after looking at the logs, I don't think the more complex models have actually reached convergence yet. I will use these new "optimums" as initial parameters again and rerun the optimisations, but I will change the rounds a little bit:
```
rounds = 4
reps = [100,80,60,50]
maxiters = [100,60,40,30]
folds = [3,2,1,1]
```
Increasing the number of iterations is going to considerably increase the run time I believe. Especially increasing it so much. I expect this will run for a couple of days at least. I will split the longest models in separate runs, and require 72hrs. Having more reps at the first steps should allow the script to explore as much as possible around the parameter space of the initial parameters, which might help in finding new local optimums. Also, looking at the parameters coming out so far, the favoured models are not really all in agreement as to whether the first event was an expansion or a bottleneck. So there is still a lot to debate.
```
sbatch --time=48:00:00 --job-name=dadi_1-1_basic nesi_dadi2D_Set1-1.sh
sbatch --time=48:00:00 --job-name=dadi_1-1_splitmig nesi_dadi2D_Set1-1.sh
sbatch --time=48:00:00 --job-name=dadi_1-2_basic nesi_dadi2D_Set1-2.sh
sbatch --time=48:00:00 --job-name=dadi_1-2_splitmig nesi_dadi2D_Set1-2.sh
sbatch --time=48:00:00 --job-name=dadi_2-1_nomig nesi_dadi2D_Set2-1.sh
sbatch --time=72:00:00 --job-name=dadi_2-1_secmig nesi_dadi2D_Set2-1.sh
sbatch --time=72:00:00 --job-name=dadi_2-1_mig nesi_dadi2D_Set2-1.sh
sbatch --time=48:00:00 --job-name=dadi_2-2_nomig nesi_dadi2D_Set2-2.sh
sbatch --time=72:00:00 --job-name=dadi_2-2_secmig nesi_dadi2D_Set2-2.sh
sbatch --time=72:00:00 --job-name=dadi_2-2_mig nesi_dadi2D_Set2-2.sh
```
Done. Let's hope I don't run out of time, these are very long to run...


#### Stats for selection outliers
###### 14.3.19
Finally back to work on this stuff after being drained of life by the kea project. So, I have done some more reading to decide exactly how to go about this, because there are many options and many programs that do different things, etc...and I have decided that the easiest and fastest way to go about this, considering my current time constraints, is to get the stats I need from VCFtools and some playing with R. Specifically, I am going to consider Fst, π, dxy and Tajima's D. These few stats aren't comprehensive, but they complement each other in ways that should allow me to interpret the patterns sensibly. At least I hope so.
I will plot these stats in comparisons between South Island and North Island, because the idea is that I am looking for specific differentiation between the two subspecies here, as opposed to a general pattern of local adaptation and isolation-by-distance that I will test with the environmental correlations and maybe EEMS.
First of all, everyone seems to agree that it is important to remove sites with very low minor allele frequency for this kind of tests. It seems to muddle up things and increase the rate of false positives. So, I will take the snps that I have already filtered for the popgen part, and add the MAF filter before proceeding here.
```bash
mkdir selection_stats
cd selection_stats/
cp ../pop_structure/maxmiss90_common_snps_HWE_LD.recode.vcf .
# checking how many snps I have in unplaced scaffolds
grep "^un" maxmiss90_common_snps_HWE_LD.recode.vcf | wc -l
979
```
I think that I will also remove the snps on unplaced scaffolds, for ease of cleanup later: since I plan on windowing all these stats and most of these scaffolds wouldn't be long enough to be included in the windows anyway.
```bash
vcftools=/usr/local/vcftools_0.1.15/bin/vcftools # since I am working in boros
awk 'BEGIN{OFS="\t";} /^un/{print $1,$2}' maxmiss90_common_snps_HWE_LD.recode.vcf > sc_to_exclude.txt
$vcftools --vcf maxmiss90_common_snps_HWE_LD.recode.vcf \
--exclude-positions sc_to_exclude.txt --maf 0.05 \
--remove-filtered-all \
--recode --recode-INFO-all \
--out filtered_snps_for_selection_tests
After filtering, kept 92 out of 92 Individuals
Outputting VCF file...
After filtering, kept 47098 out of a possible 101539 Sites
Run Time = 26.00 seconds
```
That's quite harsh, considering that only about ~1000 of those snps are excluded because of the unplaced scaffolds, I have to deduce that the rest have a tiny minor allele frequency in the overall NZ population...I wonder if the MAF is still that low in the single populations, because then they would simply be private alleles which might be quite informative in the case of local adaptation. I see this as more of a problem in the continuous pop comparison, as in the environmental correlation, than in an island comparison. I can also imagine how if the loci are actually under selection though they should be more frequent than that. Just in case, I will keep an additional set at a MAF of 0.02 and see what happens with that one.
```bash
$vcftools --vcf maxmiss90_common_snps_HWE_LD.recode.vcf \
--exclude-positions sc_to_exclude.txt --maf 0.02 \
--remove-filtered-all \
--recode --recode-INFO-all \
--out filtered_snps_for_selection_tests_maf02
After filtering, kept 92 out of 92 Individuals
Outputting VCF file...
After filtering, kept 77318 out of a possible 101539 Sites
Run Time = 40.00 seconds
```
Now, I will output the Fst stats with VCFtools on this set, before splitting it in two for the rest of the stats. BUT I need to define the population file. It seems I only need to specify the individuals present in one population and VCFtools calculates the Fst between this and all the others? --> Wrong, I need to provide two files, one for the individuals of each pop. Can get this from my ever useful population.txt file.
```bash
cp ../pop_structure/population.txt .
awk '/North/{print $1}' population.txt > NI.txt
awk '/South/{print $1}' population.txt > SI.txt
$vcftools --vcf filtered_snps_for_selection_tests.recode.vcf \
--weir-fst-pop NI.txt --weir-fst-pop SI.txt \
--out NI_vs_SI
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.025244
Weir and Cockerham weighted Fst estimate: 0.030833
After filtering, kept 47098 out of a possible 47098 Sites
Run Time = 2.00 seconds
```
Not the overall lowest Fst possible, but still pretty low, which suggests a condition very very close to panmixia. I see quite a few negative values in the output and the internet tells me that "In principle Fst scores are not impossible, as they mean that there is more variation within the population than between the two populations compared. In general, I believe it is common practice to change all the negative Fst scores to 0 and basically consider them as loci for which there is no population differentiation." from https://www.biostars.org/p/132253/ and others.
Again, this does not sound impossible to me.
I forgot to mention that I am getting all these stats on a per-site basis, because I will bin them in R. BUT, I can't do that for the Tajima's D, since it is a stat that is by default calculated over bp windows. I will use the windows from this stat to fix all the others as well, deciding which ones to exclude due to too few polymorphisms or end of a scaffold which doesn't reach the required size.
Anyway, first I will split the vcf into the two pops.
```bash
$vcftools --vcf filtered_snps_for_selection_tests.recode.vcf \
--keep NI.txt \
--remove-filtered-all \
--recode --recode-INFO-all \
--out NI_pop_for_stats
$vcftools --vcf filtered_snps_for_selection_tests.recode.vcf \
--keep SI.txt \
--remove-filtered-all \
--recode --recode-INFO-all \
--out SI_pop_for_stats
```
There are 49 individuals kept in the NI and 43 in the SI. I think the sample sizes are fairly even.
From these two files I will calculate Tajima's D, π and allele frequencies to use in R for dxy. These stats need to be population specific, that's why I am doing it this way.
```bash
$vcftools --vcf NI_pop_for_stats.recode.vcf \
--freq --out NI
$vcftools --vcf NI_pop_for_stats.recode.vcf \
--site-pi --out NI
$vcftools --vcf NI_pop_for_stats.recode.vcf \
--TajimaD 50000 --out NI
$vcftools --vcf SI_pop_for_stats.recode.vcf \
--freq --out SI
$vcftools --vcf SI_pop_for_stats.recode.vcf \
--site-pi --out SI
$vcftools --vcf SI_pop_for_stats.recode.vcf \
--TajimaD 50000 --out SI
```
Taking a quick peek at the Tajima's D values, I think I will have to get rid of quite a few bins that have no SNPs. I am still wondering if it makes sense to eliminate the MAF for this specific stat, since it is based on low frequency alleles. But most of all, probably I should increase the window size? This is were I think that being a bit less conservative on the MAF might be beneficial after all.
```bash
$vcftools --vcf filtered_snps_for_selection_tests_maf02.recode.vcf \
--weir-fst-pop NI.txt --weir-fst-pop SI.txt \
--out NI_vs_SI_maf02
Outputting Weir and Cockerham Fst estimates.
Weir and Cockerham mean Fst estimate: 0.022555
Weir and Cockerham weighted Fst estimate: 0.029324
After filtering, kept 47098 out of a possible 47098 Sites
Run Time = 2.00 seconds
```
I expected these numbers to go up rather than down with more minor alleles. But clearly there are some minor alleles in both populations still.
Ergo, probably not under selection.
```bash
$vcftools --vcf filtered_snps_for_selection_tests_maf02.recode.vcf \
--keep NI.txt \
--remove-filtered-all \
--recode --recode-INFO-all \
--out NI_pop_for_stats_maf02
$vcftools --vcf filtered_snps_for_selection_tests_maf02.recode.vcf \
--keep SI.txt \
--remove-filtered-all \
--recode --recode-INFO-all \
--out SI_pop_for_stats_maf02
# and the stats:
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--freq --out NI_maf02
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--site-pi --out NI_maf02
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--TajimaD 50000 --out NI_maf02
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--freq --out SI_maf02
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--site-pi --out SI_maf02
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--TajimaD 50000 --out SI_maf02
```
Number of SNPs per window is still pretty low in the Tajima's D calculations though, I think I do need to increase the window size or perform a sliding window or find a way of having similar numbers of SNPs in every window (changing the window size accordingly). I think I will increase the window size...I will lose resolution, but probably not too much, considering these stats are calculated on the number of segregating sites, not on the size of the window overall. I think what I could do is perform the Tajima's D on a few step sizes, then pick in R the windows with a certain number of SNPs only. Let's try that:
```bash
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--TajimaD 25000 --out NI_maf02_25kb
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--TajimaD 25000 --out SI_maf02_25kb
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--TajimaD 100000 --out NI_maf02_100kb
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--TajimaD 100000 --out SI_maf02_100kb
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--TajimaD 150000 --out NI_maf02_150kb
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--TajimaD 150000 --out SI_maf02_150kb
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--TajimaD 200000 --out NI_maf02_200kb
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--TajimaD 200000 --out SI_maf02_200kb
$vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--TajimaD 250000 --out NI_maf02_250kb
$vcftools --vcf SI_pop_for_stats_maf02.recode.vcf \
--TajimaD 250000 --out SI_maf02_250kb
```
I could probably simply count the number of SNPs per window and decide the split on that.
Just a quick check on the allele frequency files tells me that the reference allele is always output first in these files, which is very nice, because it means I have the same first allele and second allele in the two files and I can use them for calculations without troubles. I will give a quick fix to these frequency files before importing them in R, I only want numbers and not genotypes:
```bash
list="A C G T"
for i in $list
do
sed -i 's/'"${i}"'://' NI_maf02.frq
sed -i 's/'"${i}"'://' SI_maf02.frq
done
```
Perfect, the rest of this work will be in R, log in `GBS_plotting_outlier_stats.Rmd`.
Moving files on desktop, so I can work in RStudio there:
`scp -r boros:/data/denise/ModPop_analysis/selection_stats .`

###### 18.3.19
While working in R to get those bin sizes, I realised that the files I produces in VCFtools for Tajima's D would not work, because they could not possibly cover all the combinations of start and end of the bins that I have. I need a sliding window approach to actually be able to do that. So, I am installing VCF-kit to use their sliding window utility for Tajima's D.
It has quite a lot of dependencies, and needs to be installed in Python2.7.
```
source activate python2.7
pip install yahmm
pip install matplotlib
pip install VCF-kit
/Users/denisemartini/miniconda3/envs/python2.7/bin/vk tajima
```
This was still giving me errors, so I had to fix a couple of lines in one of the scripts:
```
nano /Users/denisemartini/miniconda3/envs/python2.7/lib/python2.7/site-packages/vcfkit/utils/vcf.py
import sys # added this
np.set_printoptions(threshold=sys.maxsize) # changed this
```
Now it runs, but we will need to see if it works:
```
cd /Volumes/Denise/ModPop_analysis/selection_stats/
/Users/denisemartini/miniconda3/envs/python2.7/bin/vk tajima
usage:
  vk tajima [--no-header --extra] <window-size> <step-size> <vcf>
  vk tajima [--no-header --extra] <window-size> --sliding <vcf>
```
Let's make some space and try:
```bash
rm *_maf02_*0kb.*
/Users/denisemartini/miniconda3/envs/python2.7/bin/vk tajima 50000 25000 NI_pop_for_stats_maf02.recode.vcf
KeyError: 'AC'
```
Looks like we need the AC info in the INFO field. Which is not going to work for me, since my platypus files don't have that field. Will need to go back to the B plan and calculate Tajima's D some other way, using ANGSD or something else. I will think about it, for now I will just plot the other stats.
I had a quick thought: maybe I could restrict VCFtools to operate within certain regions and give it the right window size for the Tajima's D calculations and see if it does what I want...testing out:
```bash
# put these regions in a bed file:
ps_ch_1 125000 250000
ps_ch_1 775000 900000
ps_ch_1 1150000 1275000
ps_ch_1 1275000 1400000
ps_ch_1 1625000 1750000
# then run:
vcftools --vcf NI_pop_for_stats_maf02.recode.vcf \
--bed regions.bed \
--TajimaD 125000 --out NI_maf02_125kb
```
Nope, it only uses the snps in the specified regions, but it still starts and ends at the regular chromosome positions. No good.
I think I will try PopGenome in R instead, which sounds promising.
That actually did work fine. The rest of the log on this topic is in the same outlier stats markdown document.

#### Environmental correlations
###### 19.3.19
For the environmental association tests I decided to use the R package LEA (http://membres-timc.imag.fr/Olivier.Francois/LEA/software.htm), which uses the lfmm method and has been used quite a lot. I downloaded the environmental variables from WorldClim (http://worldclim.org/version2), because they seem to be quite accurate and they include several biologically relevant variables in the bioclim section. I am also going to use the info from this tutorial (http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/environmental-correlation-analysis/#sec4) to extract the variables I need for my samples and put together the input files in R. The only thing I need to prepare then is the snp file, which they suggest should be cleaned up of maf < .05 and is easier to deal with if in .ped format. I can do this with VCFtools once again. I think I will be conservative and only remove maf < .02 for now. I will use the file I prepared for the selection stats and convert that.
```bash
cd env_correlations
VCFtools --gzvcf ../selection_stats/filtered_snps_for_selection_tests_maf02.recode.vcf.gz \
--plink --out snps_for_env_tests
```
I forgot that the plink format screws up a bit when you have chromosome names different from integers, but it does not matter, because I had used the chromosome name and position as SNP IDs, so I can retrace the snp identity from the .map file. LEA does not use that information anyway.
The rest of the work for this part is going to be in `GBS_env_association.Rmd`.

###### 22.3.19
After some initial testing, I am a bit afraid that I am getting loads of outliers. I think I might try rerunning the analysis on maf > .05.
Just to avoid some of those false positives.
```bash
cd env_correlations
VCFtools --vcf ../selection_stats/filtered_snps_for_selection_tests.recode.vcf \
--plink --out snps_for_env_tests
```
Then I am moving this to boros, together with the locations file and the rmarkdown script, which I am compiling like:
```R
rmarkdown::render("GBS_env_association.Rmd")
```
It will probably run for a full day before I see results.

###### 5.4.19
I had to think about what to do with those results, because it turns out that there's very many outliers in these tests.
To see if there is anything in there worth discussing, I think I can only do two things: first, I annotate the outliers, to see if they end up in genes and are non-synonymous substitutions or anything like that. Second, I also run an enrichment test to see if they are around genes involved in specific categories of processes.
I should already have everything I need to run SnpEff (on NeSI though) and GOwinda (on my desktop), from when I did this on the spp comparison project.
I need to copy over those reference files though and see if the names of chromosomes etc need fixing.
```bash
cd env_correlations
mkdir annotation
cp -r /Volumes/Denise/Kea-Kaka_take2/realignment/snpeff/data .
cp ../../../Kea-Kaka_take2/realignment/snpeff/gowinda/kaka_annotation.gtf .
cp ../../../Kea-Kaka_take2/realignment/snpeff/gowinda/GO_mappings .
cp ../../../Kea-Kaka_take2/realignment/snpeff/gowinda/ref_ann.txt .
```
I can already see that I will need to fix the chromosome names to fit the snp files. I also just realised that I don't have a proper vcf on which to run these tests, I converted it straight away to .ped, so I will need to fix that as well.
```bash
sed -i.tmp 's/ps_chr/ps_ch/' kaka_annotation.gtf
sed -i.tmp 's/ps_chr/ps_ch/g' data/Nmer/genes.gff
cp ../../selection_stats/filtered_snps_for_selection_tests.recode.vcf ./snps_for_env_tests.vcf
```
I will also need to go back from the snp index that lfmm gave me in the output of most significant snps, to the snp chromosome and position. It should be straightforward using the plink .map file (which is basically made for that). But I will run all this in NeSI because it is crazy how slow everything is now on my desktop.
```bash
cd ..
scp -r annotation mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/env_correlations/
scp -r sig*.txt mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/env_correlations/
scp -r *.map mahuika:/nesi/nobackup/uoo02327/denise/ModPop_analysis/env_correlations/
```
I just realised in the meantime that I made a bit of a fatal mistake in assuming that I would be able to get FDR<0.01 values from a list that came with FDR<0.05. The significant outliers in the list contain the p-values, not the corrected q-values, unfortunately. I will fix the rmarkdown to add that output file.
```bash
TAB=$'\t'
for i in $(grep -v "index" sig_ann_srad_01.txt | awk '{print $1}')
do
  awk "NR==$i" snps_for_env_tests.map | awk '{print $2}' | sed 's/\(ps_ch_[0-9A-Z]*\)_\([0-9]*\)/\1'"${TAB}"'\2/' >> pos_sig_ann_srad_01.txt
done
```
Once you substitute that file name with the other environmental variables, all files are ready.
These files can be run in GOwinda already as they are. Installing GOwinda on NeSI (it's a java jar file, so it is a matter of downloading it rather than installing it).
```bash
wget https://sourceforge.net/projects/gowinda/files/latest/download
mv download Gowinda-1.12.jar
```
So, now to run GOwinda on our set of candidates, it's enough in a quick nesi script:
```bash
gowinda=/nesi/project/uoo02327/programs/Gowinda-1.12.jar
filelist="ann_srad coldest_month ann_prec"
for f in $filelist
do
  java -Xmx12g -jar $gowinda --annotation-file kaka_annotation.gtf --gene-set-file GO_mappings \
  --snp-file snps_for_env_tests.vcf --candidate-snp-file pos_sig_${f}_01.txt \
  --output-file gowinda_${f}.txt --mode snp --gene-definition updownstream5000 --simulations 1000000 --threads 8
done
```
This command is set up to take in the file containing all snps used for environmental analysis, the kaka annotation file and the mappings of kaka genes to GO terms, and test against all these the set of candidate snps, which in this case includes the outliers from the environmental tests. Since I am pretty sure that the vast majority of these snps (coming from GBS) are not going to be in genes, I am expanding the gene definition to include 5kbp up- and down-stream of the gene regions in the annotations, to allow snps in regulatory regions to be associated with genes as well. Since I believe that the snps in this set are mostly in linkage equilibrium (I thinned the set, it comes from GBS, etc) I am using the mode "snp" for the analysis.
The documentation says that running 1M simulations should take about 30mins on 8 threads. Since I am looping this through my 3 sets of outliers, I will require ~2hrs for the run.
