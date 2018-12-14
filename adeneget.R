## I want to use adegenet to do some exploratory analysis of snp data (possibly a PCA) and maybe to extract info I need for a circos plot (specifically, I want allele frequencies per population)

library(adegenet)

## from the genomics tutorial:

# import from plink .raw files (+ .map files)
read.PLINK()

## to extract allele frequencies from genlight objects, replacing NAs
tab(x)
# so if I first subset it by population and then extract allele frequencies I should have my circos tracks ready? 
# I would get a list of allele frequencies that I could easily associate back to my snp positions

## to make PCA from genlight objects
pca1 <- glPca(x)

scatter(pca1, posi="bottomright")
title("PCA of the US influenza datann axes 1-2")

col <- funky(15)
s.class(pca1$li, pop(microbov),xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)
s.class(pca1$li,pop(microbov),xax=1,yax=3,sub="PCA 1-3",csub=2) # takes pca 1 and 3 , colours pops and adds ellipses around them
abline(v=0,h=0,col="grey", lty=2)  # the line on the 0 mark
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3) # adds the eigenvalues inset

## to make a NJ tree with genlight object:
library(ape)
tre <- nj(dist(as.matrix(flu)))
tre
plot(tre, typ="fan", cex=0.7)
title("NJ tree of the US influenza data")

## or there is the dapc
dapc1 <- dapc(x, n.pca=10, n.da=1)
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:2), col=c("red","blue"))
scatter(dapc1) # this is already the kind of plot I want

myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft") # with customization
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6)) # withouth labels, with legend


##########################################################################################

## trying out with the common biallelic snp dataset to see what I can do
test_set <- read.PLINK(file="biallelic_common_snps5_forR.raw", map.file = "biallelic_common_snps5.recode.plink.map")
# importing the .map file, there might be problems with variantID being not recognised

# it loaded up fine and in seconds (very good!) and it is a 6.8Mb genlight object
pop(test_set)
test_set@pop
# population data looking good
indNames(test_set)
# individual data also looking good

# trying some exploratory analysis
glPlot(test_set)

myFreq <- glMean(test_set)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*0.8,lwd=3)

myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)

# trying to get the populations in different genlight objects, to then get pop allele frequencies
popNames(test_set)
seppop(test_set) -> separated_pops

# then I could just get allele frequencies for each pop like this:
glMean(separated_pops$Codfish) -> Codfish_mean_allele
Codfish_mean_allele[3:25]

# with this I should be able to extract allele frequencies while also replacing missing data, but it's in a matrix -> this still keeps individuals separate, that's why it is a matrix, and where an individual has a NA it fills that with the overall allele frequency, whereas glMean only counts NAs as 0s
tab(separated_pops$Codfish) -> Codfish_allele_freq
Codfish_allele_freq[1:7,3:25]

# so then if I apply the means to these numbers, they are kind of corrected for missing data?
apply(Codfish_allele_freq, 2, mean) -> Codfish_mean_noNA_allele_freq
Codfish_mean_noNA_allele_freq[3:25]

## actually, it looks like glMean also does this correction, the difference in the results is not due to that but to the decision of including ploidy in the matter or not: is it allele frequency in individuals or allele frequency overall that I want? Allele frequency overall depends on ploidy, so in this case it's basically allele freq in individuals divided by 2
## I don't think it matters in our situation, since all individuals have the same ploidy anyway --> but, if I use the allele frequency overall the numbers go from 0 (homozygote for one allele) to 1 (homozygote for the other) and that makes more sense to me 

## anyway, if I just export these lists of frequencies out and I match them up to SNP positions, I have my circos tracks ready
# then I can just as easily extract only the positions I want, or do that straight in the circos plot --> very good!

##### time to try some pca:
pca_test <- glPca(test_set)
scatter(pca_test, posi="bottomleft")
title("PCA of the test_set data axes 1-2")

col <- funky(8)
s.class(pca_test$scores,pop(test_set),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, sub="PCA 1-2",csub=2)
s.class(pca_test$scores,pop(test_set),xax=1,yax=3,col=transp(col,.6),
        sub="PCA 1-3",csub=2, axesell=FALSE, cstar=0, cpoint=3) # takes pca 1 and 3 , colours pops and adds ellipses around them
s.class(pca_test$scores,pop(test_set),xax=2,yax=3,col=transp(col,.6),
        sub="PCA 2-3",csub=2, axesell=FALSE, cstar=0, cpoint=3)
abline(v=0,h=0,col="grey", lty=2)  # the line on the 0 mark
add.scatter.eig(pca_test$scores[1:20],nf=3,xax=1,yax=3) # adds the eigenvalues inset

## and dapc:
dapc_test <- dapc(test_set, n.pca=40, n.da=10)
scatter(dapc_test) 

library(RColorBrewer)
myCol <- c(values=brewer.pal(8, "RdYlGn"))
scatter(dapc_test, posi.da="topleft", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft", cex=1.5, clab=0, leg=TRUE) # with customization
scatter(dapc_test, scree.pca = TRUE, bg="white", pch=20, cstar=0, col=myCol, solid=.4,
        cex=3, clab=0, leg=TRUE, posi.da="topleft", posi.pca="bottomleft") # withouth labels, with legend

scatter(dapc_test,1,1, col=col, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc_test,2,2, col=col, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

scatter(dapc_test,1,2, col=col, bg="white", pch=20, cstar=0, cex=1.5, clab=0,
        scree.da=FALSE, legend=TRUE, solid=.4)

loadingplot(dapc_test$var.contr, axis=1, thres=4e-04)
loadingplot(dapc_test$var.contr, axis=2, thres=5e-04)

tab(test_set[,7771])
## --> cool! snp that is practically fixed only in Kapiti island, rare everywhere else
tab(test_set[,23884:23885])
## these two don't seem to be linked, but they are both at high frequency everywhere but Kapiti
tab(test_set[,12619])

##########################################################################################

## now using the real data set (snps from stacks aligned to the pseudochr assembly)
stacks <- read.PLINK(file="biallelic_maxmiss90_stacks_forR.raw", map.file = "biallelic_maxmiss90_stacks_output.recode.plink.map")
# importing the .map file, there might be problems with variantID being not recognised

# checking that population and individual names are fine
pop(stacks)
stacks@pop
# OK
indNames(stacks)
# OK

# just checking missing data and allele frequency distributions before starting
glPlot(stacks)

myFreq <- glMean(stacks)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)

# trying to get the populations in different genlight objects, to then get pop allele frequencies
popNames(stacks)
seppop(stacks) -> separated_pops

# getting allele frequencies for each pop for circos:
glMean(separated_pops$Codfish) -> Codfish_mean_allele
Codfish_mean_allele[3:25]
glMean(separated_pops$Fiordland) -> Fiordland_mean_allele
Fiordland_mean_allele[3:25]
glMean(separated_pops$Westland) -> Westland_mean_allele
Westland_mean_allele[3:25]
glMean(separated_pops$Nelson) -> Nelson_mean_allele
Nelson_mean_allele[3:25]
glMean(separated_pops$Kapiti) -> Kapiti_mean_allele
Kapiti_mean_allele[3:25]
glMean(separated_pops$Zealandia) -> Zealandia_mean_allele
Zealandia_mean_allele[3:25]
glMean(separated_pops$Pureora) -> Pureora_mean_allele
Pureora_mean_allele[3:25]
glMean(separated_pops$LittleBarrier) -> LittleBarrier_mean_allele
LittleBarrier_mean_allele[3:25]

## making the dapc plot:
dapc_stacks <- dapc(stacks, n.pca=40, n.da=10)
scatter(dapc_stacks) 

library(RColorBrewer)
myCol <- c("#D73027", "#F46D43", "#D9EF8B", "#1A9850", "#FEE08B", "#66BD63", "#FDAE61", "#A6D96A")
scatter(dapc_stacks, posi.da="topleft", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft", cex=1.5, clab=0, leg=TRUE) # with customization
scatter(dapc_stacks, scree.pca = TRUE, bg="white", pch=20, cstar=0, col=myCol, solid=.6,
        cex=3, clab=0, leg=TRUE, posi.da="bottomright", posi.pca="bottomleft") # withouth labels, with legend

scatter(dapc_stacks,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc_stacks,2,2, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

scatter(dapc_stacks,1,2, col=myCol, bg="white", pch=20, cstar=0, cex=1.5, clab=0,
        scree.da=FALSE, legend=TRUE, solid=.4)

loadingplot(dapc_stacks$var.contr, axis=1, thres=0.0008)
loadingplot(dapc_stacks$var.contr, axis=2, thres=0.0009)

tab(stacks[,66917])
