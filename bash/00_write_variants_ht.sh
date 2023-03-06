#!/usr/bin/env bash
#
# @description Write Hail Table of variants 
# @author Frederik Heymann Lassen (with edits by Nik Baya)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=write_variants_ht
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/write_variants_ht.log
#SBATCH --error=logs/write_variants_ht.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue
#SBATCH --array=1-23

set -x
set -o errexit
set -o nounset

for util_type in {qsub,hail}; do
  source "/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ukb_utils/bash/${util_type}_utils.sh"
done

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="python/ukb_wes_450k_gwas/00_export_csqs.py" # Can use this script with certain flags to just write a table

readonly chr=$( get_chr ${SLURM_ARRAY_TASK_ID} ) 
readonly in_dir="data/vcf"
readonly in_vcf="${in_dir}/ukb_wes_450k.qced.sites_only.chr${chr}.vcf.bgz"

readonly out_dir="data/ht"
readonly out_prefix="${out_dir}/ukb_wes_450k.qced.chr${chr}"

mkdir -p ${spark_dir} ${out_dir}

if [ ! -f "${out_prefix}.ht/_SUCCESS" ]; then
  set_up_hail
  python3 "${hail_script}" \
    --in_vcf ${in_vcf}\
    --out_prefix ${out_prefix} \
    --skip_and_write_ht_of_variants \
    && print_update "Finished writing table for chr${chr}" ${SECONDS} \
    || raise_error "Writing table for chr${chr} failed"
else
  >&2 echo "${out_prefix}.ht/_SUCCESS already exists. Skipping."
fi 

