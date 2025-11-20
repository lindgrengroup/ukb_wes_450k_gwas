#!/usr/bin/env bash
#
# Author: Nik Baya (2023-08-04)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gene_annot
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/get_imputed_v3_gene_annots.log
#SBATCH --error=logs/logs/get_imputed_v3_gene_annots.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 2
#SBATCH --mem-per-cpu=10G
#SBATCH --requeue
#SBATCH --array=1-23

set -u # Raise error if unbound variable is used

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"

# [OPTION] chrom
# Chromosome to run
if [ ${SLURM_ARRAY_TASK_ID} -eq 23 ]; then
  readonly chrom='X'
else
  readonly chrom=${SLURM_ARRAY_TASK_ID}
fi

python_script="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/python/ukb_wes_450k_gwas/get_imputed_v3_gene_annotations.py"

python3 ${python_script} --chrom ${chrom}
