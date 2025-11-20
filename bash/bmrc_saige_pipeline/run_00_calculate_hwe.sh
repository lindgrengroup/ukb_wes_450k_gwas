#!/usr/bin/env bash
#
# Author: Nik Baya (2023-03-16)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=calc_hwe
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/calculate_hwe.log
#SBATCH --error=logs/calculate_hwe.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
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


## [INPUT]
readonly bfile_dir="${WD}/data/00_set_up"
readonly bfile="${bfile_dir}/ukb_impv3.subset_to_wes_450k_qc_pass_${pop}.chr${chrom}"

## [OUTPUT]
readonly out=${bfile}

plink2 \
  --bfile ${bfile} \
  --hardy \
  --out ${bfile} \
  && gzip ${bfile}.hardy*

if [[ ${chrom} == 'X' ]]; then
  new_file=$( echo ${bfile}.hardy* | sed 's/hardy.x/hardy/g' )
  mv ${bfile}.hardy* ${new_file}
fi
