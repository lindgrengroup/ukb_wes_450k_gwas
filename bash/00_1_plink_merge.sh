#!/bin/bash

set -e # Stop job if any command fails

## OPTIONS
# pop:
# - 'allpop'
# - 'eur'
readonly pop=$1

# dataset:
# - 'pruned'
# - 'for_vr'
readonly dataset=$2 

# include_chrX: 
# - true
# - false
readonly include_chrX=true


## INPUT
get_bfile() {
	local _dataset=$1
	local _chrom=$2
	
	echo "/mnt/project/saige_pipeline/data/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.${_dataset}.chr${_chrom}"
}

## OUTPUT
merged_bfile="ukb_array.wes_450k_qc_pass_${pop}.${dataset}"


for chrom in {1..22}; do
	get_bfile ${dataset} ${chrom} >> merge_list.txt
done

if ${include_chrX}; then
	get_bfile ${dataset} "X" >> merge_list.txt
fi

plink \
	--merge-list merge_list.txt \
	--make-bed \
	--out "${merged_bfile}"

rm merge_list.txt