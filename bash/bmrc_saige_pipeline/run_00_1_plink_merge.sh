#!/usr/bin/env bash
#
# Used to merge PLINK datasets
#
# Author: Nik Baya (2023-03-09)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=merge
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/saige_merge.log
#SBATCH --error=logs/saige_merge.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 16
#SBATCH --requeue

set -e # Stop job if any command fails

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"

# [OPTION] pop
# Population to use
#   Options:
#   - 'eur': Genetically European
#   - 'allpop': Individuals from all populations who pass QC (i.e. no population filter)
readonly pop="eur"

# [OPTION] dataset
# - 'pruned' # DEPRECATED (not subset to imputed data)
# - 'subset_to_impv3.pruned'
# - 'for_vr'
readonly dataset="for_vr" 

# [OPTION] include_chrX
# - true
# - false
readonly include_chrX=true


# [INPUT]
get_bfile() {
    local _dataset=$1
    local _chrom=$2
    
    echo "${WD}/data/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.${_dataset}.chr${_chrom}"
}

# [OUTPUT]
readonly out_dir="${WD}/data/00_set_up"
readonly merged_bfile="${out_dir}/ukb_array.wes_450k_qc_pass_${pop}.${dataset}"

mkdir -p ${out_dir}

merge_list="/tmp/merge_list"
rm -f ${merge_list}

for chrom in {1..22}; do
    get_bfile ${dataset} ${chrom} >> ${merge_list}
done

if ${include_chrX}; then
    get_bfile ${dataset} "X" >> ${merge_list}
fi

plink \
  --merge-list ${merge_list} \
  --make-bed \
  --out "${merged_bfile}"

rm -f ${merge_list}
