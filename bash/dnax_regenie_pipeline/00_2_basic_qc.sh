#!/usr/bin/env bash

anc=$1 # Ancestry group

ORIGINAL_BFILE="ukb22418_b0_v2.autosomes"
PASS_QC="${ORIGINAL_BFILE}.qced.$anc"

# Subset to specific ancestry
# Necessary to avoid invariant sites in non-EUR ancestry
KEEP_WES_QC_DNAX="/mnt/project/brava/inputs/ancestry_sample_ids/qced_${anc}_sample_IDs.txt"
KEEP_WES_QC_LOCAL="$HOME/tmp-sample_ids.wes_qc_pass_$anc.txt"
awk '{ print $1,$1 }' $KEEP_WES_QC_DNAX > $KEEP_WES_QC_LOCAL

# Suggested QC filters by https://rgcgithub.github.io/regenie/recommendations/
plink2 \
  --bfile "/mnt/project/nbaya/regenie/data/genotypes/$ORIGINAL_BFILE" \
  --keep ${KEEP_WES_QC_LOCAL} \
  --maf 0.01 \
  --mac 100 \
  --geno 0.1 \
  --hwe 1e-15 0 \
  --mind 0.1 \
  --write-snplist \
  --write-samples \
  --no-id-header \
  --out $PASS_QC

