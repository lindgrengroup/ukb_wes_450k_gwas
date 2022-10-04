#!/bin/bash
#
# Create PLINK dataset with 2000 randomly selected markers, using the hard-called genotypes.
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

readonly chrom=$1
readonly array_bfile="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c${chrom}_b0_v2"
readonly out="ukb_array_for_vr_chr${chrom}"

#1. Calculate allele counts for each marker in the large PLINK file with hard called genotypes
plink2 \
  --bfile "${array_bfile}" \
  --freq counts \
  --out "${out}"

#2. Randomly extract IDs for markers falling in the two MAC categories:
# * 1,000 markers with 10 <= MAC < 20
# * 1,000 markers with MAC >= 20
cat <(
  tail -n +2 "${out}.acount" \
  | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' \
  | shuf -n 1000 ) \
<( \
  tail -n +2 "${out}.acount" \
  | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' \
  | shuf -n 1000 \
  ) > "${out}.markerid.list"


#3. Extract markers from the large PLINK file
plink2 \
  --bfile "${array_bfile}" \
  --extract "${out}.markerid.list" \
  --make-bed \
  --out "${out}"

