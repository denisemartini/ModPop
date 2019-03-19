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
kakalocs <- read.table("../env_correlations/samples_locations.txt", sep="\t", header=T)
kakalocs <- tibble(kakalocs$SampleID, kakalocs$Longitude, kakalocs$Latitude)
colnames(kakalocs) <- c("SAMPLE", "LONG", "LAT")
```

Loading the bioclimatic variables files.


```r
files <- list.files("../env_correlations/env_variables/wc2.0_2.5m_bio/", pattern='tif', full.names=TRUE)
bioclim2.5 <- stack(files)
```

Solar radiation is in a stack of its own (wd) and needs to be synthesized in a single layer if I want to have one value for the annual solar radiation, rather than twelve monthly values. I decided to go for a mean value.


```r
files_wd <- list.files("../env_correlations/env_variables/wc2.0_2.5m_srad/", pattern='tif', full.names=TRUE)
wd <- stack(files_wd)
mean(wd) -> mean_srad
```

Now I can extract the values for my locations from the variables I want. And I can plot the maps too, because it's fun. And save the result as .env files for LEA.
Annual mean temperature.


```r
mean_temp <- tibble(kakalocs$SAMPLE, extract(bioclim2.5[[1]], kakalocs[,2:3]))
plot(bioclim2.5, 1, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
write.table(mean_temp, "../env_correlations/mean_temp.env", col.names = F, row.names = F, quote = F)
```

Minimum temperature of the coldest month.


```r
coldest_month <- tibble(kakalocs$SAMPLE, extract(bioclim2.5[[6]], kakalocs[,2:3]))
plot(bioclim2.5, 6, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
write.table(coldest_month, "../env_correlations/coldest_month.env", col.names = F, row.names = F, quote = F)
```

Annual precipitation.


```r
ann_prec <- tibble(kakalocs$SAMPLE, extract(bioclim2.5[[12]], kakalocs[,2:3]))
plot(bioclim2.5, 12, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
write.table(ann_prec, "../env_correlations/ann_prec.env", col.names = F, row.names = F, quote = F)
```

Solar radiation.


```r
ann_srad <- tibble(kakalocs$SAMPLE, extract(mean_srad, kakalocs[,2:3]))
plot(mean_srad, xlim=c(160,180), ylim=c(-48,-34))
points(x=kakalocs$LONG, y=kakalocs$LAT, col="red", cex=0.8)
```

![](GBS_env_association_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
write.table(ann_srad, "../env_correlations/ann_srad.env", col.names = F, row.names = F, quote = F)
```
