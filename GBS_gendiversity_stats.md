---
title: "Working on some genetic diversity stats from the GBS data"
author: "Denise Martini"
date: "3/14/2019"
output: 
  html_document: 
    keep_md: yes
---

I will use this log to experiment on how to get the specific stats I want from the files I have, then loop this over all populations in a separate script. 


```r
library(dplyr)
```

So, first of all, I think I can get the heterozygosity measures, which at this stage I mean as average proportion of sites at which an individual of that population was heterozygous. I have a file with the per-individual homozygosity and number of sites, I just need to extract this info for all the individuals in a population and average it across them.
First, let's load the files with the individuals in a population (Kapiti here) and heterozygosity.


```r
read.table("../pop_structure/gen_diversity/Kapiti.txt", header = FALSE) -> indv
read.table("../pop_structure/gen_diversity/all_het.het", header = TRUE) -> het
```

I can very quickly take from the het file only the rows that are in common with the indv file, then attach another column that is the proportion of het on available sites (1-hom/all). The mean of this new column is the value I want.


```r
pop_het <- left_join(indv, het, by=c("V1"="INDV"))
```

```
## Warning: Column `V1`/`INDV` joining factors with different levels, coercing
## to character vector
```

```r
pop_het <- mutate(pop_het, prop_het=(1-(pop_het$O.HOM./pop_het$N_SITES)))
mean(pop_het$prop_het)
```

```
## [1] 0.1215247
```

Done. Now, a bit more complicated, I need the proportion of fixed sites and the average minor allele frequencies of the sites that are polymorphic in each population. First, loading the allele frequencies file.


```r
read.table("../pop_structure/gen_diversity/Kapiti.frq", header = TRUE, row.names = NULL) -> freq
colnames(freq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_REF", "FREQ_ALT")
nrow(freq) -> total_snps
```

Now, I need to filter out all rows with missing values in the allele frequencies, then filter out all homozygous sites. I also need the counts of the homozygous and missing sites.


```r
freq %>% filter(., is.na(FREQ_REF)) -> missing
nrow(missing) -> missing_snps
missing_snps/total_snps #proportion of missing snps
```

```
## [1] 0.01211357
```

```r
freq %>% filter(., !is.na(FREQ_REF)) -> freq
freq %>% filter(., FREQ_REF >=1 | FREQ_REF <= 0) -> hom
nrow(hom) -> hom_snps
hom_snps/total_snps #proportion of snps fixed in this population
```

```
## [1] 0.5987059
```

```r
freq %>% filter(., !(FREQ_REF >=1 | FREQ_REF <= 0)) -> freq
nrow(freq) -> het_snps
het_snps/total_snps
```

```
## [1] 0.3891805
```

Finally, the average minor allele frequency. I will need to choose the minimum value from the last two rows, then calculate the mean.


```r
apply(freq[,5:6], 1, FUN=min) -> maf
mean(maf)
```

```
## [1] 0.2118519
```

I will also output the positions of ALT hom and het sites of each population, so I can compare them and find sites that are private of one pop. I will also need fixed alleles in general for one of the things I have in mind, so outputting that as well.


```r
hom %>% filter(., FREQ_REF <= 0)  %>% select(., CHROM, POS) -> sites
bind_rows(sites, freq[,1:2]) %>% arrange(., CHROM, POS) -> sites
write.table(sites, "../pop_structure/gen_diversity/Kapiti_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(hom[,1:2], "../pop_structure/gen_diversity/Kapiti_fixed.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
```

Now I will fix the code above in a script that I can loop through my populations, I will run it and come back here for that final merge looking for private alleles. 
Okay, the script ran fine, now to look at these private alleles. I could load up all the positions, from one of the frequency files, then add the other pops as I go, like a TRUE/FALSE thing, then sum up the number of trues to find the ones that are present in only one pop.


```r
read.table("../pop_structure/gen_diversity/Kapiti.frq", header = TRUE, row.names = NULL) -> freq
freq[,1:2] -> all_sites
paste(all_sites$row.names, all_sites$CHROM, sep = "_") -> all_sites
read.table("../pop_structure/gen_diversity/Kapiti_sites.txt", header = F, row.names = NULL) -> Kapiti
paste(Kapiti$V1, Kapiti$V2, sep = "_") -> Kapiti
all_sites %in% Kapiti -> included_sites
read.table("../pop_structure/gen_diversity/LittleBarrier_sites.txt", header = F, row.names = NULL) -> LittleBarrier
paste(LittleBarrier$V1, LittleBarrier$V2, sep = "_") -> LittleBarrier
read.table("../pop_structure/gen_diversity/Pureora_sites.txt", header = F, row.names = NULL) -> Pureora
paste(Pureora$V1, Pureora$V2, sep = "_") -> Pureora
read.table("../pop_structure/gen_diversity/Zealandia_sites.txt", header = F, row.names = NULL) -> Zealandia
paste(Zealandia$V1, Zealandia$V2, sep = "_") -> Zealandia
read.table("../pop_structure/gen_diversity/Nelson_sites.txt", header = F, row.names = NULL) -> Nelson
paste(Nelson$V1, Nelson$V2, sep = "_") -> Nelson
read.table("../pop_structure/gen_diversity/Westland_sites.txt", header = F, row.names = NULL) -> Westland
paste(Westland$V1, Westland$V2, sep = "_") -> Westland
read.table("../pop_structure/gen_diversity/Fiordland_sites.txt", header = F, row.names = NULL) -> Fiordland
paste(Fiordland$V1, Fiordland$V2, sep = "_") -> Fiordland
read.table("../pop_structure/gen_diversity/Codfish_sites.txt", header = F, row.names = NULL) -> Codfish
paste(Codfish$V1, Codfish$V2, sep = "_") -> Codfish
included_sites <- tibble((all_sites %in% LittleBarrier), (all_sites %in% Pureora), (all_sites %in% Kapiti), (all_sites %in% Zealandia),
                         (all_sites %in% Nelson), (all_sites %in% Westland), (all_sites %in% Fiordland), (all_sites %in% Codfish))
```

Now that I have all that, I can try and figure out what sites (if any) are present in only one pop.


```r
rownames(included_sites) <- all_sites
tibble::rownames_to_column(included_sites) -> included_sites
filter(included_sites, included_sites[,2]==TRUE & included_sites[,3]==FALSE & included_sites[,4]==FALSE &
         included_sites[,5]==FALSE & included_sites[,6]==FALSE & included_sites[,7]==FALSE & included_sites[,8]==FALSE &
         included_sites[,9]==FALSE) %>% select(.,rowname) -> LittleBarrier_private
filter(included_sites, included_sites[,3]==TRUE & included_sites[,2]==FALSE & included_sites[,4]==FALSE &
         included_sites[,5]==FALSE & included_sites[,6]==FALSE & included_sites[,7]==FALSE & included_sites[,8]==FALSE &
         included_sites[,9]==FALSE) %>% select(.,rowname) -> Pureora_private
filter(included_sites, included_sites[,4]==TRUE & included_sites[,3]==FALSE & included_sites[,2]==FALSE &
         included_sites[,5]==FALSE & included_sites[,6]==FALSE & included_sites[,7]==FALSE & included_sites[,8]==FALSE &
         included_sites[,9]==FALSE) %>% select(.,rowname) -> Kapiti_private
filter(included_sites, included_sites[,5]==TRUE & included_sites[,3]==FALSE & included_sites[,4]==FALSE &
         included_sites[,2]==FALSE & included_sites[,6]==FALSE & included_sites[,7]==FALSE & included_sites[,8]==FALSE &
         included_sites[,9]==FALSE) %>% select(.,rowname) -> Zealandia_private
filter(included_sites, included_sites[,6]==TRUE & included_sites[,3]==FALSE & included_sites[,4]==FALSE &
         included_sites[,5]==FALSE & included_sites[,2]==FALSE & included_sites[,7]==FALSE & included_sites[,8]==FALSE &
         included_sites[,9]==FALSE) %>% select(.,rowname) -> Nelson_private
filter(included_sites, included_sites[,7]==TRUE & included_sites[,3]==FALSE & included_sites[,4]==FALSE &
         included_sites[,5]==FALSE & included_sites[,6]==FALSE & included_sites[,2]==FALSE & included_sites[,8]==FALSE &
         included_sites[,9]==FALSE) %>% select(.,rowname) -> Westland_private
filter(included_sites, included_sites[,8]==TRUE & included_sites[,3]==FALSE & included_sites[,4]==FALSE &
         included_sites[,5]==FALSE & included_sites[,6]==FALSE & included_sites[,7]==FALSE & included_sites[,2]==FALSE &
         included_sites[,9]==FALSE) %>% select(.,rowname) -> Fiordland_private
filter(included_sites, included_sites[,9]==TRUE & included_sites[,3]==FALSE & included_sites[,4]==FALSE &
         included_sites[,5]==FALSE & included_sites[,6]==FALSE & included_sites[,7]==FALSE & included_sites[,8]==FALSE &
         included_sites[,2]==FALSE) %>% select(.,rowname) -> Codfish_private
```

Let's just output the sites and numbers:

```r
nrow(LittleBarrier_private)
```

```
## [1] 3855
```

```r
write.table(LittleBarrier_private, "../pop_structure/gen_diversity/LittleBarrier_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
nrow(Pureora_private)
```

```
## [1] 1656
```

```r
write.table(Pureora_private, "../pop_structure/gen_diversity/Pureora_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
nrow(Kapiti_private)
```

```
## [1] 3542
```

```r
write.table(Kapiti_private, "../pop_structure/gen_diversity/Kapiti_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
nrow(Zealandia_private)
```

```
## [1] 3191
```

```r
write.table(Zealandia_private, "../pop_structure/gen_diversity/Zealandia_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
nrow(Nelson_private)
```

```
## [1] 1048
```

```r
write.table(Nelson_private, "../pop_structure/gen_diversity/Nelson_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
nrow(Westland_private)
```

```
## [1] 4059
```

```r
write.table(Westland_private, "../pop_structure/gen_diversity/Westland_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
nrow(Fiordland_private)
```

```
## [1] 1964
```

```r
write.table(Fiordland_private, "../pop_structure/gen_diversity/Fiordland_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
nrow(Codfish_private)
```

```
## [1] 1246
```

```r
write.table(Codfish_private, "../pop_structure/gen_diversity/Codfish_private_sites.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
```

These numbers look like a direct function of the number of individuals I had in each population. Which does make sense.
