## Modern Population study - analysis log
#### Denise Martini, 17.11.18
---

All the analysis is done in the ModPop_analysis directory in boros/nesi and moved to HCS for storage

#### What needs to be done:
- [x] Fixing the input files (genome, keyfiles, scripts, etc)
- [ ] Rerun variant callers, specifically realignments and then:
  - [ ] platypus
  - [ ] stacks
  - [ ] tassel 5
  - [ ] GATK
  - [ ] ipyrad
- [ ] Filter and merge results, vcftools and VennDiagram
- [ ] Population structure tests:
  - [ ] admixture
  - [ ] dapc in adegenet
  - [ ] tree in treemix?
  - [ ] modeling in dadi?
- [ ] Stats for selection outliers (Fst, Tajima's D, etc)
- [ ] Environmental correlations
- [ ] Inbreeding?

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

Finally, moving the script to run Tassel5 into the repo and fixing it
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

I had forgotten an important detail: Tassel runs BWA at some point for the alignment to the reference .fasta genome, and it assumes that this is indexed already. In this case it wasn't, so it failed to locate the index and did not run the alignment, with obvious consequences for the rest of the pipeline. So:
```bash
cd ..
/usr/local/bwa-0.7.17/bwa index pseudochromosomes.fasta
cd tassel
bash tassel5-GBS2.sh 2>&1 | tee tassel.log
```
