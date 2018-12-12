#!/bin/bash -e

## 13.12.18, Denise
# fixed from a previous script, in order to run with different random seeds as well as different Ks
# set variables
admixture=/usr/local/admixture_1.3.0/admixture
seedlist="76 11 24 161 358"
Klist="1 2 3 4 5 6 7 8 9 10"

for s in $seedlist
do
  for K in $Klist
  do
    $admixture -s $s --cv maxmiss90_common_snps.ped $K -j8 | tee admix_${s}_${K}.log
  done
done
