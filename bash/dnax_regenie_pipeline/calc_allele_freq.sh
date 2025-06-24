#!/usr/bin/env bash

anc=$1
PHENOCOL=$2
chrom=$3

bfile="ukb_wes_450k.qced.chr${chrom}"
bfile_path="/mnt/project/Barney/wes/sample_filtered/${bfile}"

# Use DNAnexus:nbaya/regenie/notebooks/get_sample_list.ipynb to create ID list
ids_w_defined_pheno="/mnt/project/nbaya/regenie/data/allele_freq/ukb_wes_450k.qced.${anc}.${PHENOCOL}.iids.txt"
keep=tmp-keep.txt
awk '{ print 0, $1 }' ${ids_w_defined_pheno} > $keep

out="${bfile}.${anc}.${PHENOCOL}"
# if [ $( dx ls /nbaya/regenie/data/allele_freq/${out}.frq | wc -l 2> /dev/null ) -eq 0 ]; then
plink \
  --bfile ${bfile_path} \
  --keep $keep \
  --freq \
  --out ${out}

rm *log tmp* *nosex
