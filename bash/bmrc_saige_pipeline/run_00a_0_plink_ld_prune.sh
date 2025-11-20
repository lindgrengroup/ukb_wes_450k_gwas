#!/usr/bin/env bash
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#
# Author: Nik Baya (2023-03-08)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=ld_prune
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/saige_ld_prune.log
#SBATCH --error=logs/saige_ld_prune.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 4
#SBATCH --requeue
#SBATCH --array=1-23

set -e # Stop job if any command fails

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


## [INPUT]
readonly bed="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr${chrom}_v2.bed"
readonly bim="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr${chrom}_v2.bim"
readonly fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
readonly samples_w_superpop="${WD}/data/00_set_up/ukb_wes_450k.qced.subset_to_impv3.sample_list_w_superpops.tsv"

## [OUTPUT]
readonly out_dir="${WD}/data/00_set_up"
readonly out="${out_dir}/ukb_array.wes_450k_qc_pass_${pop}.subset_to_impv3.pruned.chr${chrom}"

mkdir -p ${out_dir}

if [[ "${pop}" == "eur" ]]; then
  # Subset to a genetic ancestry
  # LD prune on European subset of samples passing WES 450k QC
  plink2 \
    --bed "${bed}" \
    --bim "${bim}" \
    --fam "${fam}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --indep-pairwise 50 5 0.05 \
    --out "${out}"

  # Extract set of pruned variants and export to bfile
  plink2 \
    --bed "${bed}" \
    --bim "${bim}" \
    --fam "${fam}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --extract "${out}.prune.in" \
    --make-bed \
    --out "${out}"
fi
