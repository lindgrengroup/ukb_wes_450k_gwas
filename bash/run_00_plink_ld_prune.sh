#!/bin/bash

chrom=21

dx build 00_plink_ld_prune --overwrite

dx run 00_plink_ld_prune \
  -igenotype_array_bfile="wes_450k:/Bulk/Genotype Results/Genotype calls/ukb22418_c#CHROM#_b0_v2" \
  -isamples_w_superpop="wes_450k:/saige_pipeline/data/00_set_up/ukb_wes_450k.qced.sample_list_w_superpops.tsv" \
  --destination "/tmp/" \
  -y --brief \
  --priority="low" \
  --watch

