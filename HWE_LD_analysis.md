---
title: "Investigation of HWE and LD per population"
author: "Denise Martini"
date: "12/12/2018"
output: 
  html_document: 
    keep_md: yes
---



In this script I am gathering the per population statistics that VCFtools output for Hardy-Weinberg Equilibrium and Linkage Disequilibrium, to identify loci that are out of HWE and/or in LD in more than 3 populations. 


```r
library(dplyr)
```

#### HWE

First, trialling with one population only, I want to import the file, select only the columns I need (position and p-value). The file contains the calculated p-value at ALL sites, so I also need to filter only loci with a p-value<0.05.


```r
kapiti_hwe <- read.table("../vcf_filtering/Kapiti.hwe", header = TRUE) 
kapiti_hwe <- select(kapiti_hwe, CHR, POS, P_HWE)
kapiti_hwe <- filter(kapiti_hwe, P_HWE<0.05)
head(kapiti_hwe)
```

```
##       CHR      POS        P_HWE
## 1 ps_ch_1   835317 0.0107622000
## 2 ps_ch_1   882637 0.0303030300
## 3 ps_ch_1  4163917 0.0476190500
## 4 ps_ch_1  5823700 0.0052173910
## 5 ps_ch_1 11657092 0.0007876469
## 6 ps_ch_1 15960165 0.0476190500
```

That was easy enough, so writing it out in a loop that does the above filtering for each population file and adds the output to the same table.


```r
hwe_table <- NULL
for (pop in c(list.files(path = "../vcf_filtering/", pattern = ".hwe", full.names=TRUE))) {
  pop_hwe <- read.table(pop, header = TRUE)
  pop_hwe <- select(pop_hwe, CHR, POS, P_HWE)
  pop_hwe <- filter(pop_hwe, P_HWE<0.05)
  
  if(is.null(hwe_table)) {
    hwe_table <- pop_hwe
  } else {
    hwe_table <- full_join(hwe_table, pop_hwe, by=c("CHR", "POS"))
  }
}
```

I then just need to filter only the loci that are out of HWE (have a p-value rather than NA) in 3 or more populations. 


```r
hwe_table %>% mutate(., popnum=apply(hwe_table[3:(dim(hwe_table)[2])], 
                                     MARGIN = 1, function(x){ sum(!is.na(x)) })) -> hwe_table
out_of_hwe <- filter(hwe_table, popnum>=3)
out_of_hwe <- arrange(out_of_hwe, CHR, POS)
head(out_of_hwe)
```

```
##       CHR      POS    P_HWE.x     P_HWE.y    P_HWE.x.x    P_HWE.y.y
## 1 ps_ch_1   835317         NA 0.031674210 0.0107622000 0.0044420510
## 2 ps_ch_1 11657092 0.03729604 0.006906406 0.0007876469 0.0004939601
## 3 ps_ch_1 19403268 0.03729604 0.002903186 0.0007876469 0.0004939601
## 4 ps_ch_1 28145481         NA 0.009287926           NA           NA
## 5 ps_ch_1 28178990 0.03729604 0.010530650 0.0056912100 0.0043463300
## 6 ps_ch_1 49429438         NA          NA           NA 0.0274605500
##   P_HWE.x.x.x P_HWE.y.y.y P_HWE.x.x.x.x P_HWE.y.y.y.y popnum
## 1          NA          NA  1.962579e-02   0.004524887      5
## 2          NA 0.006906406  2.112463e-04   0.001856402      7
## 3          NA 0.006906406  5.616715e-05   0.001856402      7
## 4          NA 0.046439630  3.893215e-02            NA      3
## 5          NA          NA  2.610718e-03            NA      5
## 6          NA 0.010530650  2.610718e-03   0.014287640      4
```

```r
dim(out_of_hwe)[1] #to check the number of loci found, count the number of rows
```

```
## [1] 237
```

Looking good, writing out the IDs of the filtered loci in a format that VCFtools accepts.


```r
write.table(out_of_hwe[,1:2], file = "../vcf_filtering/loci_outofHWE.txt", sep = "_", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

#### LD

The input file format here is a bit different, it already only contains significant sites (with an rsq>0.5). There is a column that has the number of individuals that were used for the LD calculation at those two positions (all individuals that had a genotype for both loci). At many positions only few individuals were used for the calculation and I don't think that's right. So, I am filtering out values for which less than 6 individuals were used in the comparisons. After that, I can get rid of that column. Again, trialling first with one population.


```r
kapiti_ld <- read.table("../vcf_filtering/Kapiti.geno.ld", header = TRUE)
kapiti_ld <- filter(kapiti_ld, N_INDV>=6)
kapiti_ld <- select(kapiti_ld, c(CHR, POS1, POS2, R.2))
head(kapiti_ld)
```

```
##       CHR  POS1   POS2      R.2
## 1 ps_ch_1 10540  29295 0.606061
## 2 ps_ch_1 10540  60927 0.606061
## 3 ps_ch_1 10540 117771 0.675000
## 4 ps_ch_1 10540 158890 1.000000
## 5 ps_ch_1 29295  60927 1.000000
## 6 ps_ch_1 29295 158890 0.606061
```

Everything looks okay, I can loop through the populations and add all results to the same table. 


```r
ld_table <- NULL
for (pop in c(list.files(path = "../vcf_filtering/", pattern = ".geno.ld", full.names=TRUE))) {
  pop_ld <- read.table(pop, header = TRUE)
  pop_ld <- filter(pop_ld, N_INDV>=6)
  pop_ld <- select(pop_ld, c(CHR, POS1, POS2, R.2))
  
  if(is.null(ld_table)) {
    ld_table <- pop_ld
  } else {
    ld_table <- full_join(ld_table, pop_ld, by=c("CHR", "POS1", "POS2"))
  }
}
```

Again, as before, I am selecting the loci that are in LD in 3 or more populations.


```r
ld_table %>% mutate(., popnum=apply(ld_table[4:(dim(ld_table)[2])], 
                                    MARGIN = 1, function(x){ sum(!is.na(x)) })) -> ld_table
in_ld <- filter(ld_table, popnum>=3)
in_ld <- arrange(in_ld, CHR, POS1)
head(in_ld)
```

```
##       CHR    POS1    POS2  R.2.x R.2.y  R.2.x.x  R.2.y.y R.2.x.x.x
## 1 ps_ch_1 1162958 1163085     NA     1 0.729167       NA        NA
## 2 ps_ch_1 3982198 4015622     NA    NA       NA 0.611111        NA
## 3 ps_ch_1 5423167 5423945 1.0000     1       NA 1.000000        NA
## 4 ps_ch_1 5519131 5520065 0.5625     1       NA 1.000000        NA
## 5 ps_ch_1 6214675 6239865     NA     1       NA 1.000000        NA
## 6 ps_ch_1 7650802 7733884     NA     1       NA 1.000000        NA
##   R.2.y.y.y R.2.x.x.x.x R.2.y.y.y.y popnum
## 1        NA    0.611111          NA      3
## 2        NA    1.000000    1.000000      3
## 3        NA          NA    1.000000      4
## 4        NA    0.786332    1.000000      5
## 5        NA          NA    0.650000      3
## 6        NA          NA    0.900901      3
```

```r
dim(in_ld)[1] #to check the number of loci found, count the number of rows
```

```
## [1] 2504
```

Now, there is a bit more tweaking to be done on these loci, because I don't need to throw them all out, I only need to discard one of two loci in LD between them. First, I rework the dataset so that the loci are identified by ID and I put all the IDs in a separate dataframe.


```r
in_ld %>% mutate_at(., .vars = "POS1", funs(paste(in_ld$CHR,in_ld$POS1, sep = "_"))) %>% 
  mutate_at(., .vars = "POS2", funs(paste(in_ld$CHR,in_ld$POS2, sep = "_"))) %>% 
  select(., CHR, POS1, POS2)-> pos_in_ld
pos1 <- select(pos_in_ld, POS1)
pos2 <- select(pos_in_ld, POS2) 
colnames(pos1) <- "POS"
colnames(pos2) <- "POS"
pos <- union_all(pos1, pos2)
```

Then, I verify that some IDs in this dataset are duplicated (loci in LD with more than one other locus). I identify all the IDs present more than once in the dataset and I store them in a separate list. 


```r
pos <- arrange(pos, POS)
pos %>% n_distinct() # number of distinct loci that are in LD with at least another locus
```

```
## [1] 3617
```

```r
pos %>% add_count(.,POS) %>% filter(.,n>1) -> dups
dups <- distinct(dups)
dups <- as.list(dups$POS)
length(dups) # number of distinct loci that are in LD with more than one other locus
```

```
## [1] 817
```

I set some rules to pick one locus from a LD pair: 

- if they are both unique loci, randomly sample one of the two

- if one of the two is a duplicate, pick the other

- if they are both duplicates, put aside for later


```r
keep <- NULL
trouble <- NULL
for (i in 1:nrow(pos_in_ld)) {
  if(!(pos_in_ld$POS1[i] %in% dups) && !(pos_in_ld$POS2[i] %in% dups)) {
    keep[i] <- sample(pos_in_ld[i,2:3],1) #sample at random
  } else if(!(pos_in_ld$POS1[i] %in% dups) && (pos_in_ld$POS2[i] %in% dups)) {
    keep[i] <- pos_in_ld$POS1[i] #pick locus1
  } else if((pos_in_ld$POS1[i] %in% dups) && !(pos_in_ld$POS2[i] %in% dups)) {
    keep[i] <- pos_in_ld$POS2[i] #pick locus2
  } else { 
    keep[i] <- NA #space filler
    trouble <- c(trouble, rownames(pos_in_ld[i,])) #record the line number of duplicates pair for later
  }
}
keep <- keep[!is.na(keep)]
length(keep)
```

```
## [1] 1612
```

The pairs that get filtered out because both loci are duplicated belong to groups of loci that are all in LD between each other, so picking one from each pair would be problematic. But, picking one locus from each group should be fine. 


```r
trouble <- as.numeric(trouble)
keeptoo <- NULL
for (i in 1:length(trouble)) {
  if (i==1) {
    keeptoo <- c(keeptoo, pos_in_ld$POS1[(trouble[i])])
  } else if (!(trouble[i]==(trouble[i-1]+1)))
    keeptoo <- c(keeptoo, pos_in_ld$POS1[(trouble[i])])
}
n_distinct(keeptoo) == length(keeptoo) #sanity check, there should be no duplicates in this set
```

```
## [1] TRUE
```

```r
keep <- c(keep, keeptoo)
length(keep)
```

```
## [1] 1823
```

Now that I know which loci I can keep, I can discard all the others.

```r
pos %>% filter(., !(POS %in% keep)) -> discard
discard <- distinct(discard)
discard <- as.list(discard$POS)
length(discard) # number of loci that will need to be discarded
```

```
## [1] 1794
```

```r
n_distinct(pos) == (length(keep)+length(discard))  # sanity check
```

```
## [1] TRUE
```

Looking good, writing out the IDs of the loci to filter in a format that VCFtools accepts. Also writing out the IDs of the retained loci, in case they are needed later on.


```r
write.table(discard, file = "../vcf_filtering/loci_inLD.txt", sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(keep, file = "../vcf_filtering/loci_inLD_kept.txt", sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
```


