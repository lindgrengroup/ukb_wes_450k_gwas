#!/usr/bin/env bash
#
# @description Write Hail Table of variants 
# @author Frederik Heymann Lassen (with edits by Nik Baya)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=concat_vcfs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/concat_vcfs.log
#SBATCH --error=logs/concat_vcfs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

set -x
set -o errexit
set -o nounset

for util_type in {qsub,hail}; do
  source "/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ukb_utils/bash/${util_type}_utils.sh"
done

module load BCFtools/1.12-GCC-10.3.0

chrom=$1
max_parts=$2

file_list=/tmp/bcftools_file_list.txt

vcf_prefix="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/vcf/ukb_wes_450k.qced.sites_only.chr"
echo "${vcf_prefix}${chrom}-1of${max_parts}.vcf.bgz" > ${file_list}

for i in `seq 2 ${max_parts}`; do
  echo "${vcf_prefix}${chrom}-${i}of${max_parts}.vcf.bgz" >> ${file_list}
done

bcftools concat \
 --file-list ${file_list} \
 | bgzip > "${vcf_prefix}${chrom}.vcf.bgz"
