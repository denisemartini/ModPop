## Modern Population study - analysis log
#### Denise Martini, 17.11.18
---

All the analysis is done in the ModPop_analysis directory in boros/nesi and moved to HCS for storage

#### What needs to be done:
- [ ] Fixing the input files (genome, keyfiles, scripts, etc)
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
