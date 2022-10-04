#!/bin/bash

readonly out="ukb_array_pruned"

get_bfile() {
	local chrom=$1

	echo "/mnt/project/saige_pipeline/data/00_sparse_grm/ukb_array_pruned_chr${chrom}"
}

for chrom in {1..22}; do
	get_bfile ${chrom} >> merge_list.txt
done
get_bfile "X" >> merge_list.txt

plink \
	--merge-list merge_list.txt \
	--make-bed \
	--out ${out}

rm merge_list.txt