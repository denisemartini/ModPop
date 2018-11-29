## Modern Population study - analysis log
#### Denise Martini, 17.11.18
---

All the analysis is done in the ModPop_analysis directory in boros/nesi and moved to HCS for storage

#### What needs to be done:
- [x] Fixing the input files (genome, keyfiles, scripts, etc)
- [ ] Rerun variant callers, specifically realignments and then:
  - [x] platypus
  - [x] stacks
  - [x] tassel 5
  - [x] GATK
  - [ ] ipyrad
- [ ] Filter and merge results, vcftools and VennDiagram
- [ ] Population structure tests:
  - [ ] admixture
  - [ ] dapc in adegenet
  - [ ] tree in treemix
  - [ ] modeling in dadi
- [ ] Stats for selection outliers (Fst, Tajima's D, etc)
- [ ] Environmental correlations
- [ ] Inbreeding tests

###### _Technical note_
Most of the first part of the analysis has already been tested before, and scripts/tips are available from the Kaka_GBS directory. Parts yet to test from that part are ipyrad, treemix, dadi.  
One thing that would be worth investigating is if I can run most of this in NeSI this time around, and if it would be faster that way. From a quick check, it should be easy enough to run the alignments, stacks, ipyrad, GATK in NeSI, since they are already installed in there. For the demultiplexing, sabre would need to be added. Platypus and Tassel would also need to be installed, and while Platypus was easy enough (but also fast enough that it would probably be easy to run in boros), Tassel was a nightmare for Hugh the first time around, so I might just run it in boros either way. It was slow, so should be started as quickly as possible. Admixture would probably be an easy install, most of the rest runs in R. Need to investigate dadi.

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
