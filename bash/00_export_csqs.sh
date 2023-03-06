#!/usr/bin/env bash
#
# @description Get variant csqs and MAF/MAC count by combined tables
# @author Frederik Heymann Lassen (with edits by Nik Baya)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=export_csqs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/export_csqs.log
#SBATCH --error=logs/export_csqs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-23
#SBATCH --requeue

set -o errexit
set -o nounset

for util_type in {bash,qsub,hail}; do
    source "/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ukb_utils/bash/${util_type}_utils.sh"
  done

readonly spark_dir="data/tmp/spark_dir"
readonly hail_script="python/00_export_csqs.py"
readonly rscript="rscripts/00_export_csqs.R"

readonly chr=$( get_chr ${SGE_TASK_ID} ) 
readonly by="worst_csq_by_gene_canonical"
readonly in_dir="data/vcf"
readonly in_vcf="${in_dir}/ukb_wes_450k.qced.chr${chr}.vcf.bgz"

readonly out_dir="data/vep/${by}"
readonly out_prefix="${out_dir}/ukb_wes_450k.qced.chr${chr}"
readonly out_saige="${out_prefix}.saige"

mkdir -p ${spark_dir} ${out_dir}

if [ ! -f "${out_prefix}.tsv.gz" ]; then
  set_up_hail
  set_up_pythonpath_legacy  
  python3 "${hail_script}" \
    --in_vcf ${in_vcf}\
    --out_prefix ${out_prefix} \
    --by "${by}" \
    --by_explode \
    --out_type "tsv" \
    && print_update "Finished exporting csqs chr${chr}" ${SECONDS} \
    || raise_error "Exporting csqs for chr${chr} failed"
else
  >&2 echo "${out_prefix}.tsv.gz already exists. Skipping.."
fi 


# Generate SAIGE-GENE+ Group file consequence
# annotations (SAIGE version > 0.99.2)
module purge
set_up_rpy
Rscript ${rscript} \
  --input_path "${out_prefix}.tsv.gz" \
  --output_path "${out_saige}" \
  --delimiter " "

