#!/bin/bash

set -e # Stop job if any command fails

dataset=$1 # Options: pruned, for_vr
include_chrX=true

merged_bfile="ukb_array.wes_450k_qc_pass_eur.${dataset}"

get_bfile() {
	local chrom=$1
	local dataset=$2

	echo "/mnt/project/saige_pipeline/data/00_set_up/ukb_array.wes_450k_qc_pass_eur.${dataset}.chr${chrom}"
}

for chrom in {1..22}; do
	get_bfile ${chrom} ${dataset} >> merge_list.txt
done
if ${include_chrX}; then
	get_bfile "X" >> merge_list.txt
fi

plink \
	--merge-list merge_list.txt \
	--make-bed \
	--out "${merged_bfile}"

rm merge_list.txt