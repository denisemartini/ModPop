suppressMessages(library(dplyr, quietly=T))

#usage: Rscript GBS_gendiverity_stats.R <popfile> <hetfile>
#eg: Rscript GBS_gendiverity_stats.R poplist.txt all_het.het

#get arguments from the command line
args<-commandArgs(TRUE)
stopifnot(length(args)==2)
if(length(args)==0){
  print("No arguments specified; cannot proceed")
} else if(length(args)==1){
  print("Either poplist or heterozygosity file are missing; cannot proceed")
} else if(length(args)==2){
  poplist<-as.character(args[1])
  hetfile<-as.character(args[2])
} else {
  print("Too many arguments given from the command line")
}

print(paste0("Reading input files"))
# getting input files
read.table(poplist, header = F, sep = "\t") -> pops
read.table(hetfile, header = T, sep = "\t") -> het

# start variables to fill
meanhet <- c()
missing_snps <- c()
hom_snps <- c()
het_snps <- c()
meanmaf <- c()

# loop through pops
for (i in seq(1,to=(nrow(pops)),by=1)) {
  popname=pops[i,1]
  print(paste0("Working on pop ", popname))
  indfile=(paste0(popname,".txt")) # read in file with individuals for that pop
  read.table(indfile, header = FALSE) -> indv
  frqfile=(paste0(popname,".frq")) # read in frequency file for that pop
  read.table(frqfile, header = TRUE, row.names = NULL) -> freq
  colnames(freq) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_REF", "FREQ_ALT")
  
  # mean heterozygosity
  pop_het <- left_join(indv, het, by=c("V1"="INDV"))
  pop_het <- mutate(pop_het, prop_het=(1-(pop_het$O.HOM./pop_het$N_SITES)))
  meanhet <- c(meanhet, mean(pop_het$prop_het))
  
  # filter out freq file
  freq %>% filter(., is.na(FREQ_REF)) -> missing
  c(missing_snps, nrow(missing)) -> missing_snps
  freq %>% filter(., !is.na(FREQ_REF)) -> freq
  freq %>% filter(., FREQ_REF >=1 | FREQ_REF <= 0) -> hom
  c(hom_snps, nrow(hom)) -> hom_snps
  freq %>% filter(., !(FREQ_REF >=1 | FREQ_REF <= 0)) -> freq
  c(het_snps, nrow(freq)) -> het_snps
  
  # calculate mean maf
  apply(freq[,5:6], 1, FUN=min) -> maf
  c(meanmaf, mean(maf)) -> meanmaf
  
  # output files for this pop
  hom %>% filter(., FREQ_REF <= 0)  %>% select(., CHROM, POS) -> sites
  bind_rows(sites, freq[,1:2]) %>% arrange(., CHROM, POS) -> sites
  outsites=paste0(popname,"_sites.txt")
  outfixed=paste0(popname,"_fixed.txt")
  write.table(sites, outsites, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(hom[,1:2], outfixed, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

# putting together all other variables in one output
nrow(freq) -> total_snps
outstats <- tibble(pops[,1],meanhet,missing_snps,hom_snps,het_snps,meanmaf)
write.table(outstats, "gendiversity_stats.txt", col.names = c("POP","MEAN_HET","MISSING","HOM","HET","MEAN_MAF"), 
            row.names = FALSE, sep = "\t", quote = FALSE)