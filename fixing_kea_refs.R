library(dplyr)

setwd("/Volumes/Denise/ModPop_analysis/pop_structure/dadi/")

# importing the snps for both kaka and kea, these are very simple and clean vcf files so there's no need for special imports
read.table("noheader_kaka_GBS_snps.vcf", header = T, sep = "\t", comment.char = "") -> full_kaka_snps 
read.table("noheader_kea_GBS_snps.vcf", header = T, sep = "\t", comment.char = "") -> full_kea_snps 

# getting only the portions that I need to compare
kaka_snps <- select(full_kaka_snps, X.CHROM, POS, REF)
kea_snps <- select(full_kea_snps, X.CHROM, POS, REF)

# finding the sites at which kea has a different ref than kaka
setdiff(kea_snps, kaka_snps) -> problem_rows
# the sites where all is fine
intersect(kea_snps, kaka_snps) -> okay_rows

# now that I have checked the references, I only need the positions
problem_rows <- select(problem_rows, X.CHROM, POS)

# I can take the correct ref and alt alleles from kaka
with_alt_kaka_snps <- select(full_kaka_snps, X.CHROM, POS, REF, ALT)
inner_join(problem_rows, with_alt_kaka_snps) -> refs_to_fix

# select good and bad sites from the full kea set
inner_join(problem_rows, full_kea_snps) -> kea_to_fix
inner_join(okay_rows, full_kea_snps) -> kea_fine_rows

# at this stage I visually checked that the error was always the same: 
# ref and alt allele got somehow inverted in the variant calling in kea
# fix the bad sites, substituting the ref and alt alleles with the correct ones from kaka
kea_to_fix <- c(kea_to_fix[,1:3], refs_to_fix[,3:4], kea_to_fix[6:10])
as.data.frame(kea_to_fix) -> kea_to_fix

# fix the genotypes, inverting 1/1 and 0/0, because I technically inverted ref and alt alleles
kea_to_fix %>% filter(., kea=="0/0") %>% mutate(., kea = replace(kea, kea=="0/0", "1/1")) -> kea_hom_0
kea_to_fix %>% filter(., kea=="1/1") %>% mutate(., kea = replace(kea, kea=="1/1", "0/0")) -> kea_hom_1
# the rest (heterozygous) don't need fixing
kea_to_fix %>% filter(., !(kea=="1/1" | kea=="0/0")) -> kea_other

# put all the fixed bits together again, with the rest, and sort it
bind_rows(kea_hom_0, kea_hom_1, kea_other, kea_fine_rows) -> all_kea_fixed
arrange(all_kea_fixed, X.CHROM, POS) -> all_kea_fixed

# write out a new vcf file
write.table(all_kea_fixed, "fixed_kea_GBS_snps", sep = '\t', quote = F, row.names = F, col.names = F)
