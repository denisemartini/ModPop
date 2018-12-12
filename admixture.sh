
# set variables
admixture=/usr/local/admixture_1.3.0/admixture
Klist="1 2 3 4 5 6 7 8 9 10"

for K in $Klist
do

$admixture --cv ../biallelic_common_snps.bed $K -j8 | tee log${K}.out

done
