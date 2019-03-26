---
title: "Environmental association analysis on GBS data"
author: "Denise Martini"
date: "19 March 2019"
output: 
  html_document: 
    keep_md: yes
---

I want to use the package LEA to run some environmental association tests on my GBS data. First, I need to setup the environmental data though. I downloaded the variables I want to use from the WorldClim dataset, and I have a file with the coordinates (longitude and latitude) for my kaka individuals (or rather, their populations).
Now, I needed to install a few packages, but that went quite smoothly, they are here for reference. 


```r
install.packages("dismo")
install.packages("maptools")
install.packages("rgdal")
```

Now, loading the packages that I actually need for this first part.


```r
library(dismo)
library(maptools)
library(dplyr)
```

Loading and fixing my locations file, because I realised the order of the coordinates was wrong. These packages like to have longitude before latitude.


```r
kakalocs <- read.table("samples_locations.txt", sep="\t", header=T)
kakalocs <- tibble(kakalocs$SampleID, kakalocs$Longitude, kakalocs$Latitude)
colnames(kakalocs) <- c("SAMPLE", "LONG", "LAT")
```

Loading the bioclimatic variables files.


```r
files <- list.files("env_variables/wc2.0_2.5m_bio/", pattern='tif', full.names=TRUE)
bioclim2.5 <- stack(files)
```

Solar radiation is in a stack of its own (wd) and needs to be synthesized in a single layer if I want to have one value for the annual solar radiation, rather than twelve monthly values. I decided to go for a mean value.


```r
files_wd <- list.files("env_variables/wc2.0_2.5m_srad/", pattern='tif', full.names=TRUE)
wd <- stack(files_wd)
mean(wd) -> mean_srad
```

Now I can extract the values for my locations from the variables I want. And I can plot the maps too, because it's fun. And save the result as .env files for LEA.
Annual mean temperature.


```r
mean_temp <- extract(bioclim2.5[[1]], kakalocs[,2:3])
plot(bioclim2.5, 1, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
write.table(mean_temp, "mean_temp.env", col.names = F, row.names = F, quote = F)
```

Minimum temperature of the coldest month.


```r
coldest_month <- extract(bioclim2.5[[6]], kakalocs[,2:3])
plot(bioclim2.5, 6, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
write.table(coldest_month, "coldest_month.env", col.names = F, row.names = F, quote = F)
```

Annual precipitation.


```r
ann_prec <- extract(bioclim2.5[[12]], kakalocs[,2:3])
plot(bioclim2.5, 12, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
write.table(ann_prec, "ann_prec.env", col.names = F, row.names = F, quote = F)
```

Solar radiation.


```r
ann_srad <- extract(mean_srad, kakalocs[,2:3])
plot(mean_srad, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
write.table(ann_srad, "ann_srad.env", col.names = F, row.names = F, quote = F)
```

Now, on to use LEA for the correlations, real quick.
Nope, not on this version of R. -.-
So, sicne I can't update R on my desktop, because the OS is too out of date, I will finish writing up this script and run it in boros instead, where I have my nice R 3.5 installed and I already checked that LEA works fine.


```r
install.packages("devtools")
devtools::install_github("bcm-uga/LEA")
```

Loading packages and snp files, that I need to quickly convert in other formats.


```r
library(LEA)
library(RColorBrewer)
snpfile = ped2lfmm("snps_for_env_tests.ped")
```

```
## 
## 	- number of detected individuals:	92
## 	- number of detected loci:		47098
```

```r
genotype = lfmm2geno("snps_for_env_tests.lfmm")
```

```
## 
## 	- number of detected individuals:	92
## 	- number of detected loci:		47098
```

Testing the population structure (this is a similar thing to admixture), just to check that the signal is the same as usual.


```r
obj.snmf = snmf(genotype, K = 1:8, entropy = T, ploidy = 2, project="new", CPU = 8)
```

```
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1631232403
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1631232403
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 1  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          1
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K1/run1/snps_for_env_tests_r1.1.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K1/run1/snps_for_env_tests_r1.1.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1631232403
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 
## Least-square error: 1634756.886086
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K1/run1/snps_for_env_tests_r1.1.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K1/run1/snps_for_env_tests_r1.1.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      1
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K1/run1/snps_for_env_tests_r1.1.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K1/run1/snps_for_env_tests_r1.1.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.638844
## Cross-Entropy (masked data):	 0.685418
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] "*************************************"
## [1] "* sNMF K = 2  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          2
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K2/run1/snps_for_env_tests_r1.2.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K2/run1/snps_for_env_tests_r1.2.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1631232403
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[========================]
## Number of iterations: 63
## 
## Least-square error: 1586271.762042
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K2/run1/snps_for_env_tests_r1.2.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K2/run1/snps_for_env_tests_r1.2.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      2
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K2/run1/snps_for_env_tests_r1.2.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K2/run1/snps_for_env_tests_r1.2.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.60861
## Cross-Entropy (masked data):	 0.672893
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] "*************************************"
## [1] "* sNMF K = 3  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          3
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K3/run1/snps_for_env_tests_r1.3.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K3/run1/snps_for_env_tests_r1.3.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  2015970894227
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[======================================]
## Number of iterations: 102
## 
## Least-square error: 1550570.915972
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K3/run1/snps_for_env_tests_r1.3.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K3/run1/snps_for_env_tests_r1.3.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      3
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K3/run1/snps_for_env_tests_r1.3.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K3/run1/snps_for_env_tests_r1.3.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.585926
## Cross-Entropy (masked data):	 0.669246
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1631232403
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[==============]
## Number of iterations: 37
## 
## Least-square error: 1519456.423131
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565312
## Cross-Entropy (masked data):	 0.669583
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] "*************************************"
## [1] "* sNMF K = 5  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          5
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K5/run1/snps_for_env_tests_r1.5.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K5/run1/snps_for_env_tests_r1.5.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1631232403
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[===========================================================================]
## Number of iterations: 200
## 
## Least-square error: 1495640.791947
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K5/run1/snps_for_env_tests_r1.5.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K5/run1/snps_for_env_tests_r1.5.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      5
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K5/run1/snps_for_env_tests_r1.5.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K5/run1/snps_for_env_tests_r1.5.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.555088
## Cross-Entropy (masked data):	 0.671201
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] "*************************************"
## [1] "* sNMF K = 6  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          6
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K6/run1/snps_for_env_tests_r1.6.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K6/run1/snps_for_env_tests_r1.6.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1631232403
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[===========================================================================]
## Number of iterations: 200
## 
## Least-square error: 1474292.449114
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K6/run1/snps_for_env_tests_r1.6.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K6/run1/snps_for_env_tests_r1.6.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      6
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K6/run1/snps_for_env_tests_r1.6.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K6/run1/snps_for_env_tests_r1.6.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.545043
## Cross-Entropy (masked data):	 0.68361
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] "*************************************"
## [1] "* sNMF K = 7  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          7
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K7/run1/snps_for_env_tests_r1.7.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K7/run1/snps_for_env_tests_r1.7.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  4645498201635987859
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[===================]
## Number of iterations: 50
## 
## Least-square error: 1451133.543218
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K7/run1/snps_for_env_tests_r1.7.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K7/run1/snps_for_env_tests_r1.7.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      7
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K7/run1/snps_for_env_tests_r1.7.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K7/run1/snps_for_env_tests_r1.7.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.533334
## Cross-Entropy (masked data):	 0.687654
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] "*************************************"
## [1] "* sNMF K = 8  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          8
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K8/run1/snps_for_env_tests_r1.8.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K8/run1/snps_for_env_tests_r1.8.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1631232403
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[===========================================================================]
## Number of iterations: 200
## 
## Least-square error: 1429124.533319
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K8/run1/snps_for_env_tests_r1.8.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K8/run1/snps_for_env_tests_r1.8.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      8
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K8/run1/snps_for_env_tests_r1.8.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K8/run1/snps_for_env_tests_r1.8.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.523996
## Cross-Entropy (masked data):	 0.689634
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
```

```r
plot(obj.snmf)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
barplot(t(Q(obj.snmf, K = 4)), col = c(values=brewer.pal(4, "RdYlGn")))
```

![](GBS_env_association_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

Everything seems to be as usual, so I am running the imputation for missing data, to check that everything is fine again and because they specify that it is recommended before the environmental associations.


```r
obj.snmf = snmf("snps_for_env_tests.lfmm", K = 4,
                       entropy = TRUE, repetitions = 10,
                       project = "new", CPU = 8)
```

```
## [1] "*********************************************************"
## [1] " Conversion from the lfmm format to the geno format"
## [1] "*********************************************************"
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1810620828
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1810620828
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 1      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  4611686020238008732
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[===================]
## Number of iterations: 51
## 
## Least-square error: 1519492.560806
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run1/snps_for_env_tests_r1.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565442
## Cross-Entropy (masked data):	 0.668543
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1891696672
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1891696672
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 2      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run2/snps_for_env_tests_r2.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run2/snps_for_env_tests_r2.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  10481631264
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[==============]
## Number of iterations: 37
## 
## Least-square error: 1520142.563584
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run2/snps_for_env_tests_r2.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run2/snps_for_env_tests_r2.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run2/snps_for_env_tests_r2.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run2/snps_for_env_tests_r2.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565492
## Cross-Entropy (masked data):	 0.665577
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 140216186
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      140216186
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 3      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run3/snps_for_env_tests_r3.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run3/snps_for_env_tests_r3.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  140216186
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[==========================]
## Number of iterations: 69
## 
## Least-square error: 1520861.234609
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run3/snps_for_env_tests_r3.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run3/snps_for_env_tests_r3.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run3/snps_for_env_tests_r3.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run3/snps_for_env_tests_r3.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565365
## Cross-Entropy (masked data):	 0.669374
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1115907141
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1115907141
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 4      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run4/snps_for_env_tests_r4.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run4/snps_for_env_tests_r4.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1115907141
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[=====================]
## Number of iterations: 55
## 
## Least-square error: 1520269.208519
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run4/snps_for_env_tests_r4.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run4/snps_for_env_tests_r4.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run4/snps_for_env_tests_r4.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run4/snps_for_env_tests_r4.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565496
## Cross-Entropy (masked data):	 0.665883
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1506716242
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1506716242
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 5      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run5/snps_for_env_tests_r5.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run5/snps_for_env_tests_r5.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  13826050857534138962
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[==================]
## Number of iterations: 47
## 
## Least-square error: 1520528.079404
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run5/snps_for_env_tests_r5.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run5/snps_for_env_tests_r5.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run5/snps_for_env_tests_r5.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run5/snps_for_env_tests_r5.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565473
## Cross-Entropy (masked data):	 0.667598
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1805863061
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1805863061
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 6      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run6/snps_for_env_tests_r6.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run6/snps_for_env_tests_r6.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1805863061
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[===========]
## Number of iterations: 29
## 
## Least-square error: 1520704.502273
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run6/snps_for_env_tests_r6.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run6/snps_for_env_tests_r6.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run6/snps_for_env_tests_r6.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run6/snps_for_env_tests_r6.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565556
## Cross-Entropy (masked data):	 0.666079
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1006612994
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1006612994
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 7      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run7/snps_for_env_tests_r7.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run7/snps_for_env_tests_r7.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1006612994
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[========================]
## Number of iterations: 64
## 
## Least-square error: 1519855.501336
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run7/snps_for_env_tests_r7.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run7/snps_for_env_tests_r7.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run7/snps_for_env_tests_r7.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run7/snps_for_env_tests_r7.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565218
## Cross-Entropy (masked data):	 0.67049
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1168608476
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1168608476
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 8      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run8/snps_for_env_tests_r8.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run8/snps_for_env_tests_r8.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  5463575772
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[=================]
## Number of iterations: 45
## 
## Least-square error: 1519931.328634
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run8/snps_for_env_tests_r8.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run8/snps_for_env_tests_r8.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run8/snps_for_env_tests_r8.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run8/snps_for_env_tests_r8.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565503
## Cross-Entropy (masked data):	 0.66763
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1593196739
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1593196739
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 9      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run9/snps_for_env_tests_r9.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run9/snps_for_env_tests_r9.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1593196739
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[=============]
## Number of iterations: 35
## 
## Least-square error: 1520486.968151
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run9/snps_for_env_tests_r9.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run9/snps_for_env_tests_r9.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run9/snps_for_env_tests_r9.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run9/snps_for_env_tests_r9.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565406
## Cross-Entropy (masked data):	 0.67352
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
## 
## [1] 1392197177
## [1] "*************************************"
## [1] "*          create.dataset            *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)                 92
##         -L (number of loci)                        47098
##         -s (seed random init)                      1392197177
##         -r (percentage of masked data)             0.05
##         -x (genotype file in .geno format)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -o (output file in .geno format)           /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
## 
##  Write genotype file with masked data, /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## [1] "*************************************"
## [1] "* sNMF K = 4  repetition 10      *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)             92
##         -L (number of loci)                    47098
##         -K (number of ancestral pops)          4
##         -x (input file)                        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         -q (individual admixture file)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run10/snps_for_env_tests_r10.4.Q
##         -g (ancestral frequencies file)        /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run10/snps_for_env_tests_r10.4.G
##         -i (number max of iterations)          200
##         -a (regularization parameter)          10
##         -s (seed random init)                  1392197177
##         -e (tolerance error)                   1E-05
##         -p (number of processes)               8
##         - diploid
## 
## Read genotype file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno:		OK.
## 
## 
## Main algorithm:
## 	[                                                                           ]
## 	[====================]
## Number of iterations: 54
## 
## Least-square error: 1518678.549973
## Write individual ancestry coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run10/snps_for_env_tests_r10.4.Q:		OK.
## Write ancestral allele frequency coefficient file /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run10/snps_for_env_tests_r10.4.G:	OK.
## 
## [1] "*************************************"
## [1] "*    cross-entropy estimation       *"
## [1] "*************************************"
## summary of the options:
## 
##         -n (number of individuals)         92
##         -L (number of loci)                47098
##         -K (number of ancestral pops)      4
##         -x (genotype file)                 /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.geno
##         -q (individual admixture)          /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run10/snps_for_env_tests_r10.4.Q
##         -g (ancestral frequencies)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/K4/run10/snps_for_env_tests_r10.4.G
##         -i (with masked genotypes)         /data/denise/ModPop_analysis/env_correlations_maf05/snps_for_env_tests.snmf/masked/snps_for_env_tests_I.geno
##         - diploid
## 
## Cross-Entropy (all data):	 0.565411
## Cross-Entropy (masked data):	 0.667679
## The project is saved into :
##  snps_for_env_tests.snmfProject 
## 
## To load the project, use:
##  project = load.snmfProject("snps_for_env_tests.snmfProject")
## 
## To remove the project, use:
##  remove.snmfProject("snps_for_env_tests.snmfProject")
```

```r
# select the run with the lowest cross-entropy value
best = which.min(cross.entropy(obj.snmf, K = 4))
barplot(t(Q(obj.snmf, K = 4, run = best)), col = c(values=brewer.pal(4, "RdYlGn")))
```

![](GBS_env_association_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
# Impute the missing genotypes
impute(obj.snmf, "snps_for_env_tests.lfmm",
       method = 'mode', K = 4, run = best)
```

```
## Missing genotype imputation for K = 4 
## Missing genotype imputation for run = 2 
## Results are written in the file:  snps_for_env_tests.lfmm_imputed.lfmm
```

```r
file.rename("snps_for_env_tests.lfmm_imputed.lfmm", "snps_for_env_tests_imputed.lfmm")
```

```
## [1] TRUE
```

Now, time to finally run those association tests on this dataset. Starting with the mean temp variable.
In the documentation they recommend running the test with different Ks around the value identified with snmf, with multiple repetitions for each K (repetitions=5-10), a decent number of cycles (e.g. -i 6000) and half of that as burning cycles (-b 3000). So, here we go.


```r
mean_temp = lfmm("snps_for_env_tests_imputed.lfmm", "mean_temp.env", 
                K = 2:5, repetitions = 5, iterations = 6000, burnin = 3000, project="new", CPU = 8)
```

```
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K2/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5416616676
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.96	 DIC: 4320091.459 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K3/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5416616676
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.798	 DIC: 4333014.336 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K4/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5416616676
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.442	 DIC: 4333071.865 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K5/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5416616676
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.579	 DIC: 4331238.619 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K2/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           6001046476
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.47	 DIC: 4322360.893 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K3/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           6001046476
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.062	 DIC: 4333027.335 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K4/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           6001046476
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.392	 DIC: 4333072.667 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K5/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           6001046476
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.073	 DIC: 4268643.266 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K2/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5108667628
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.451	 DIC: 4321875.504 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K3/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5108667628
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.291	 DIC: 4333027.043 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K4/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5108667628
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.735	 DIC: 4333066.467 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K5/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5108667628
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.014	 DIC: 4270749.754 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K2/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5413433777
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.406	 DIC: 4322283.224 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K3/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5413433777
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.867	 DIC: 4333049.024 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K4/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5413433777
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.591	 DIC: 4333076.393 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K5/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5413433777
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.519	 DIC: 4330660.381 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K2/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5254869450
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.813	 DIC: 4317131.853 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K3/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5254869450
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.051	 DIC: 4333054.081 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K4/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5254869450
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.806	 DIC: 4333052.243 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_mean_temp.lfmm/K5/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5254869450
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              mean_temp.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	mean_temp.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.765	 DIC: 4330842.218 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_mean_temp.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_mean_temp.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_mean_temp.lfmmProject")
```

This takes a VERY LONG time.
Then I need to adjust the p-values combining the repetitions from each run, and choose the appropriate run based on the shape of the p-value histogram, which "should be flat, with a peak close to zero".


```r
p_2 = lfmm.pvalues(mean_temp, K = 2)
p_3 = lfmm.pvalues(mean_temp, K = 3)
p_4 = lfmm.pvalues(mean_temp, K = 4)
p_5 = lfmm.pvalues(mean_temp, K = 5)
par(mfrow = c(2,2))
hist(p_2$pvalues, col = "#D73027")
hist(p_3$pvalues, col = "#D9EF8B")
hist(p_4$pvalues, col = "#1A9850")
hist(p_5$pvalues, col = "#FDAE61")
```

![](GBS_env_association_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

Once confirmed that K=4 is as usual our best guess, we can go on and plot the p-values of that run and see if we have any outliers. There isn't really that much of a difference between the different Ks.


```r
plot(-log10(p_4$pvalues), pch = 19, col = "#1A9850", cex = .7)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Looks like there might be outliers here.
Then let's apply a multiple testing correction (FDR with Benjamini-Hochberg method) and extract the candidates.


```r
L = length(p_4$pvalues)
q = .05
# return a list of candidates with expected FDR .05
# Benjamini-Hochberg's algorithm:
w = which(sort(p_4$pvalues) < q * (1:L) / L)
candidates = tibble(order(p_4$pvalues)[w], sort(p_4$pvalues)[w])
dim(candidates)
```

```
## [1] 1415    2
```

```r
write.table(candidates, "sig_mean_temp.txt", sep = "\t", row.names = F,
            col.names = c("index", "p-value"), quote = F)
```

Running it all for the other variables as well. Annual precipitation.


```r
ann_prec = lfmm("snps_for_env_tests_imputed.lfmm", "ann_prec.env", 
                 K = 2:5, repetitions = 5, iterations = 6000, burnin = 3000, project="new", CPU = 8)
```

```
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K2/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           14972282714319816
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.676	 DIC: 4333045.211 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K3/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           4611686019378047944
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.075	 DIC: 4333074.37 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K4/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           950660040
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.836	 DIC: 4330378.84 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K5/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           950660040
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.281	 DIC: 4325594.35 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K2/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           906113965
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.723	 DIC: 4333043.58 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K3/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           4611686019333501869
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.917	 DIC: 4333040.605 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K4/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           906113965
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.075	 DIC: 4330261.658 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K5/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           906113965
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.02	 DIC: 4328361.178 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K2/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1335988917
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.501	 DIC: 4333042.767 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K3/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1335988917
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.525	 DIC: 4333023.105 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K4/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1335988917
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.885	 DIC: 4331168.939 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K5/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1335988917
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.073	 DIC: 4325216.48 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K2/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           148916695
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.47	 DIC: 4333033.313 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K3/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           148916695
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.839	 DIC: 4333054.272 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K4/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           148916695
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.666	 DIC: 4330443.693 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K5/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           4443883991
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.709	 DIC: 4325236.121 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K2/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1649942288
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.392	 DIC: 4333045.732 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K3/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5944909584
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED: 4333048.3	 DIC: 4333063.039 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K4/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5944909584
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.187	 DIC: 4331436.422 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_prec.lfmm/K5/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1649942288
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_prec.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_prec.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.991	 DIC: 4325605.207 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_prec.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_prec.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_prec.lfmmProject")
```

```r
p_2 = lfmm.pvalues(ann_prec, K = 2)
p_3 = lfmm.pvalues(ann_prec, K = 3)
p_4 = lfmm.pvalues(ann_prec, K = 4)
p_5 = lfmm.pvalues(ann_prec, K = 5)
par(mfrow = c(2,2))
hist(p_2$pvalues, col = "#D73027")
hist(p_3$pvalues, col = "#D9EF8B")
hist(p_4$pvalues, col = "#1A9850")
hist(p_5$pvalues, col = "#FDAE61")
```

![](GBS_env_association_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

Checking outliers with K=4 and writing out candidates.


```r
plot(-log10(p_4$pvalues), pch = 19, col = "#1A9850", cex = .7)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

```r
L = length(p_4$pvalues)
q = .05
# return a list of candidates with expected FDR .05
# Benjamini-Hochberg's algorithm:
w = which(sort(p_4$pvalues) < q * (1:L) / L)
candidates = tibble(order(p_4$pvalues)[w], sort(p_4$pvalues)[w])
dim(candidates)
```

```
## [1] 1476    2
```

```r
write.table(candidates, "sig_ann_prec.txt", sep = "\t", row.names = F,
            col.names = c("index", "p-value"), quote = F)
```

Then minimum temperature of the coldest month.


```r
coldest_month = lfmm("snps_for_env_tests_imputed.lfmm", "coldest_month.env", 
                 K = 2:5, repetitions = 5, iterations = 6000, burnin = 3000, project="new", CPU = 8)
```

```
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K2/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           833984085
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333051.391	 DIC: 4329540.497 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K3/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           833984085
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.408	 DIC: 4333043.072 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K4/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           833984085
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.142	 DIC: 4333072.16 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K5/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           833984085
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.012	 DIC: 4329815.113 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K2/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           776023149
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333051.916	 DIC: 4327973.205 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K3/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           5070990445
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.521	 DIC: 4333080.526 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K4/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           776023149
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.588	 DIC: 4333070.866 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K5/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           776023149
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.326	 DIC: 4328816.416 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K2/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1768692020
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333052.077	 DIC: 4327000.289 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K3/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1768692020
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.457	 DIC: 4333047.012 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K4/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           6063659316
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.982	 DIC: 4333078.089 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K5/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1768692020
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.003	 DIC: 4330339.075 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K2/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1437897652
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333051.733	 DIC: 4330560.397 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K3/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1437897652
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.886	 DIC: 4333078.781 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K4/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1437897652
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.932	 DIC: 4333070.191 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K5/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1437897652
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.357	 DIC: 4319836.832 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K2/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           35989878
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333051.729	 DIC: 4331364.596 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K3/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           35989878
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.774	 DIC: 4333054.482 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K4/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           35989878
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.83	 DIC: 4333062.712 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_coldest_month.lfmm/K5/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           35989878
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              coldest_month.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	coldest_month.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED: 4333049.4	 DIC: 4325347.271 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_coldest_month.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_coldest_month.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_coldest_month.lfmmProject")
```

```r
p_2 = lfmm.pvalues(coldest_month, K = 2)
p_3 = lfmm.pvalues(coldest_month, K = 3)
p_4 = lfmm.pvalues(coldest_month, K = 4)
p_5 = lfmm.pvalues(coldest_month, K = 5)
par(mfrow = c(2,2))
hist(p_2$pvalues, col = "#D73027")
hist(p_3$pvalues, col = "#D9EF8B")
hist(p_4$pvalues, col = "#1A9850")
hist(p_5$pvalues, col = "#FDAE61")
```

![](GBS_env_association_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

Checking outliers with K=4 and writing out candidates.


```r
plot(-log10(p_4$pvalues), pch = 19, col = "#1A9850", cex = .7)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

```r
L = length(p_4$pvalues)
q = .05
# return a list of candidates with expected FDR .05
# Benjamini-Hochberg's algorithm:
w = which(sort(p_4$pvalues) < q * (1:L) / L)
candidates = tibble(order(p_4$pvalues)[w], sort(p_4$pvalues)[w])
dim(candidates)
```

```
## [1] 1070    2
```

```r
write.table(candidates, "sig_coldest_month.txt", sep = "\t", row.names = F,
            col.names = c("index", "p-value"), quote = F)
```

Finally mean annual solar radiation.


```r
ann_srad = lfmm("snps_for_env_tests_imputed.lfmm", "ann_srad.env", 
                 K = 2:5, repetitions = 5, iterations = 6000, burnin = 3000, project="new", CPU = 8)
```

```
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K2/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           511197898
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333051.316	 DIC: 4330367.756 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run1/snps_for_env_tests_imputed_r1_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K3/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           511197898
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.827	 DIC: 4333082.954 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run1/snps_for_env_tests_imputed_r1_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K4/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           511197898
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.108	 DIC: 4333033.379 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run1/snps_for_env_tests_imputed_r1_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 1  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K5/run1/snps_for_env_tests_imputed_r1
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           4607275912159314634
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.305	 DIC: 4330088.051 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run1/snps_for_env_tests_imputed_r1_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K2/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1522008265
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.728	 DIC: 4331939.616 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run2/snps_for_env_tests_imputed_r2_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K3/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1522008265
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.366	 DIC:  4333055.7 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run2/snps_for_env_tests_imputed_r2_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K4/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1522008265
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.574	 DIC: 4333061.099 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run2/snps_for_env_tests_imputed_r2_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 2  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K5/run2/snps_for_env_tests_imputed_r2
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1522008265
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.821	 DIC: 4269734.774 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run2/snps_for_env_tests_imputed_r2_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K2/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1554088801
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.102	 DIC: 4331659.247 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run3/snps_for_env_tests_imputed_r3_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K3/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1554088801
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.087	 DIC: 4332981.249 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run3/snps_for_env_tests_imputed_r3_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K4/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1554088801
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.682	 DIC: 4333022.118 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run3/snps_for_env_tests_imputed_r3_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 3  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K5/run3/snps_for_env_tests_imputed_r3
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           4608342229027750753
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.905	 DIC: 4234239.895 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run3/snps_for_env_tests_imputed_r3_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K2/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1469693115
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.882	 DIC: 4332338.524 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run4/snps_for_env_tests_imputed_r4_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K3/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1469693115
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.363	 DIC: 4333073.035 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run4/snps_for_env_tests_imputed_r4_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K4/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1469693115
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333047.522	 DIC: 4333064.763 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run4/snps_for_env_tests_imputed_r4_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 4  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K5/run4/snps_for_env_tests_imputed_r4
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1469693115
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.269	 DIC: 4262628.341 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run4/snps_for_env_tests_imputed_r4_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 2  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   2
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K2/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           16103676205919623
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.222	 DIC: 4333021.126 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K2/run5/snps_for_env_tests_imputed_r5_s1.2.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 3  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   3
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K3/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1272245639
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333049.633	 DIC: 4333052.23 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K3/run5/snps_for_env_tests_imputed_r5_s1.3.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 4  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   4
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K4/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           1272245639
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333048.148	 DIC: 4332961.403 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K4/run5/snps_for_env_tests_imputed_r5_s1.4.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## [1] "********************************"
## [1] "* K = 5  repetition 5  d = 1   *"
## [1] "********************************"
## Summary of the options:
## 
##         -n (number of individuals)      92
##         -L (number of loci)             47098
##         -K (number of latent factors)   5
##         -o (output file)                snps_for_env_tests_imputed_ann_srad.lfmm/K5/run5/snps_for_env_tests_imputed_r5
##         -i (number of iterations)       6000
##         -b (burnin)                     3000
##         -s (seed random init)           14418124880539015
##         -p (number of processes (CPU))  8
##         -x (genotype file)              snps_for_env_tests_imputed.lfmm
##         -v (variable file)              ann_srad.env
##         -D (number of covariables)      1
##         -d (the dth covariable)         1
## 
## Read variable file:
##  	ann_srad.env		OK.
## 
## Read genotype file:
##  	snps_for_env_tests_imputed.lfmm		OK.
## 
## <<<<
## 	 Analyse for variable 1
## 
## 		Start of the Gibbs Sampler algorithm.
## 
## 	[                                                                           ]
## 	[===========================================================================]
## 
## 		End of the Gibbs Sampler algorithm.
## 
## 	ED:4333050.332	 DIC: 4237676.973 
## 
## 	The statistics for the run are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.dic.
## 
## 	The zscores for variable 1 are registered in:
##  		snps_for_env_tests_imputed_ann_srad.lfmm/K5/run5/snps_for_env_tests_imputed_r5_s1.5.zscore.
## 	The columns are: zscores, -log10(p-values), p-values.
## 
## 	-------------------------
## 	The execution for variable 1 worked without error.
## >>>>
## 
## The project is saved into :
##  snps_for_env_tests_imputed_ann_srad.lfmmProject 
## 
## To load the project, use:
##  project = load.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
## 
## To remove the project, use:
##  remove.lfmmProject("snps_for_env_tests_imputed_ann_srad.lfmmProject")
```

```r
p_2 = lfmm.pvalues(ann_srad, K = 2)
p_3 = lfmm.pvalues(ann_srad, K = 3)
p_4 = lfmm.pvalues(ann_srad, K = 4)
p_5 = lfmm.pvalues(ann_srad, K = 5)
par(mfrow = c(2,2))
hist(p_2$pvalues, col = "#D73027")
hist(p_3$pvalues, col = "#D9EF8B")
hist(p_4$pvalues, col = "#1A9850")
hist(p_5$pvalues, col = "#FDAE61")
```

![](GBS_env_association_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

Checking outliers with K=4 and writing out candidates.


```r
plot(-log10(p_4$pvalues), pch = 19, col = "#1A9850", cex = .7)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

```r
L = length(p_4$pvalues)
q = .05
# return a list of candidates with expected FDR .05
# Benjamini-Hochberg's algorithm:
w = which(sort(p_4$pvalues) < q * (1:L) / L)
candidates = tibble(order(p_4$pvalues)[w], sort(p_4$pvalues)[w])
dim(candidates)
```

```
## [1] 1226    2
```

```r
write.table(candidates, "sig_ann_srad.txt", sep = "\t", row.names = F,
            col.names = c("index", "p-value"), quote = F)
```
