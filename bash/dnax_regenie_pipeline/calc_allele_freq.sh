#!/usr/bin/env bash

anc=$1
pheno_group=$2
pheno_col=$3
chrom=$4

bfile="ukb_wes_450k.qced.chr${chrom}"
bfile_path="/mnt/project/Barney/wes/sample_filtered/${bfile}"

# Use DNAnexus:nbaya/regenie/notebooks/get_sample_list.ipynb to create ID list
IIDS_DIR="/mnt/project/nbaya/regenie/data/allele_freq"
#ids_w_defined_pheno="${IIDS_DIR}/ukb_wes_450k.qced.${anc}.${PHENOCOL}.iids.txt"
ids_w_defined_pheno="${IIDS_DIR}/ukb_wes_450k.qced.$anc.${pheno_group}.${pheno_col}.iids.txt"
keep=tmp-keep.txt
awk '{ print 0, $1 }' ${ids_w_defined_pheno} > $keep

out="${bfile}.${anc}.${pheno_group}.${pheno_col}"
# if [ $( dx ls /nbaya/regenie/data/allele_freq/${out}.frq | wc -l 2> /dev/null ) -eq 0 ]; then
plink \
  --bfile ${bfile_path} \
  --keep $keep \
  --freq \
  --out ${out}

rm *log tmp* *nosex
