#!/usr/bin/env bash
#
# Create PLINK dataset with 2000 randomly selected markers, using the hard-called genotypes.
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#
# Author: Nik Baya (2023-03-08)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=plink_for_vr
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/saige_plink_for_vr.log
#SBATCH --error=logs/saige_plink_for_vr.errors.log
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
#   - 'allpop': Individuals from all populations who pass QC (i.e. no population filter)
readonly pop="eur"

# [OPTION] chrom
# Chromosome to run
if [ ${SLURM_ARRAY_TASK_ID} -eq 23 ]; then
  readonly chrom='X'
else
  readonly chrom=${SLURM_ARRAY_TASK_ID}
fi


# [INPUT]
readonly bed="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr${chrom}_v2.bed"
readonly bim="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr${chrom}_v2.bim"
readonly fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
readonly samples_w_superpop="${WD}/data/00_set_up/ukb_wes_450k.qced.sample_list_w_superpops.tsv"


# [OUTPUT]
readonly out_dir="${WD}/data/00_set_up"
readonly out="${out_dir}/ukb_array.wes_450k_qc_pass_${pop}.for_vr.chr${chrom}"

mkdir -p ${out_dir}

#1. Calculate allele counts for each marker in the large PLINK file with hard called genotypes
if [[ "${pop}" == "allpop" ]]; then 
  plink2 \
    --bed "${bed}" \
    --bim "${bim}" \
    --fam "${fam}" \
    --freq counts \
    --out "${out}"

elif [[ "${pop}" == "eur" ]]; then 
  plink2 \
    --bed "${bed}" \
    --bim "${bim}" \
    --fam "${fam}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --freq counts \
    --out "${out}"
fi

#2. Randomly extract IDs for markers falling in the two MAC categories:
# * 1,000 markers with 10 <= MAC < 20
# * 1,000 markers with MAC >= 20
cat <(
  tail -n +2 "${out}.acount" \
  | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' \
  | shuf -n 1000 ) \
<( \
  tail -n +2 "${out}.acount" \
  | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' \
  | shuf -n 1000 \
  ) > "${out}.markerid.list"


#3. Extract markers from the large PLINK file
if [[ "${pop}" == "allpop" ]]; then 
  plink2 \
    --bed "${bed}" \
    --bim "${bim}" \
    --fam "${fam}" \
    --extract "${out}.markerid.list" \
    --make-bed \
    --out "${out}"

elif [[ "${pop}" == "eur" ]]; then 
  # Make sure to still subset to Europeans
  plink2 \
    --bed "${bed}" \
    --bim "${bim}" \
    --fam "${fam}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --extract "${out}.markerid.list" \
    --make-bed \
    --out "${out}"
fi
