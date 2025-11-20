#!/bin/bash
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

set -e # Stop job if any command fails

readonly chrom=$1
readonly bfile="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c${chrom}_b0_v2"
readonly out="ukb_array.wes_450k_qc_pass_allpop.pruned.chr${chrom}"

# LD prune on samples passing WES 450k QC (all populations, not just subset to European)
plink \
  --bfile "${bfile}" \
  --indep-pairwise 50 5 0.05 \
  --out "${out}"

# Extract set of pruned variants and export to bfile
plink \
  --bfile "${bfile}" \
  --extract "${out}.prune.in" \
  --make-bed \
  --out "${out}"