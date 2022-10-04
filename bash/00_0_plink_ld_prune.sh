#!/bin/bash
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

readonly chrom=$1
readonly bfile="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c${chrom}_b0_v2"
readonly out="ukb_array_pruned_chr${chrom}"

plink \
  --bfile "${bfile}" \
  --indep-pairwise 50 5 0.05 \
  --out ${out}

plink \
  --bfile "${bfile}" \
  --extract "${out}.prune.in" \
  --make-bed \
  --out ${out}
