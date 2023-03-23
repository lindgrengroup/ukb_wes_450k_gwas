#!/usr/bin/env bash
#
# Author: Nik Baya (2023-03-14)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=subset_bgen
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/subset_bgen.log
#SBATCH --error=logs/subset_bgen.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 10
#SBATCH --mem-per-cpu=10G
#SBATCH --requeue
#SBATCH --array=1-23

set -u # Raise error if unbound variable is used

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"

# [OPTION] pop
# Population to use 
#   Options: 
#   - 'eur': Genetically European
readonly pop="eur"

# [OPTION] chrom
# Chromosome to run
if [ ${SLURM_ARRAY_TASK_ID} -eq 23 ]; then
  readonly chrom='X'
else
  readonly chrom=${SLURM_ARRAY_TASK_ID}
fi

# [OPTION] output_format
# Which format to output after subsetting samples
# - "plink"
# - "bgen"
readonly output_format="bgen"

if [ "${output_format}" == "plink" ]; then
  output_format_flag="--make-bed"
elif [ "${output_format}" == "bgen" ]; then
  output_format_flag="--export bgen-1.3"
fi

## [INPUT]
readonly bgen_prefix="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_outlier_phens/data/regenie/genotypes/ukb_imputed_v3/minmaf_0.001/info0.8/chr${chrom}"
readonly samples_w_superpop="${WD}/data/00_set_up/ukb_wes_450k.qced.subset_to_impv3.sample_list_w_superpops.tsv"

## [OUTPUT]
readonly out_dir="${WD}/data/00_set_up"
readonly tmp_dir="${WD}/data/tmp"
readonly out="${out_dir}/ukb_impv3.subset_to_wes_450k_qc_pass_${pop}.chr${chrom}"

mkdir -p ${out_dir} ${tmp_dir}

keep=${tmp_dir}/tmp-keep-chr${chrom}

if [[ "${pop}" == "eur" ]]; then
  awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" > ${keep}
  
  plink2 \
    --bgen ${bgen_prefix}.bgen \
    --sample ${bgen_prefix}.sample \
    --keep ${keep} \
    ${output_format_flag} \
    --out ${out} \
  && rm ${keep} \
  || rm ${keep}
fi

if [ "${output_format}" == "bgen" ]; then
  module load BGEN/1.1.6-GCCcore-7.3.0
  bgenix \
    -g "${out}.bgen" \
    -index 
fi
