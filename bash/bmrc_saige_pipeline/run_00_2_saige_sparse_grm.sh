#!/usr/bin/env bash
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#
# Author: Nik Baya (2023-03-09)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=sparse_grm
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/saige_sparse_grm.log
#SBATCH --error=logs/saige_sparse_grm.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task=48
#SBATCH --requeue

set -e # Stop job if any command fails

readonly WD="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"

# [OPTION] pop
# Population to use
#   Options:
#   - 'eur': Genetically European
#   - 'allpop': Individuals from all populations who pass QC (i.e. no population filter)
readonly pop="eur"

# [INPUT]
readonly bfile="ukb_array.wes_450k_qc_pass_${pop}.pruned"

# [OUTPUT]
readonly out_dir="${WD}/data/00_set_up"
readonly out="long-$(( SLURM_CPUS_PER_TASK ))-ukb_array.wes_450k_qc_pass_${pop}.pruned"

mkdir -p ${out_dir}

singularity exec \
  --bind ${WD}/data/00_set_up:/mnt \
  /apps/singularity/saige_1.1.6.3.sif \
  createSparseGRM.R  \
    --plinkFile="/mnt/${bfile}" \
    --nThreads=$(( SLURM_CPUS_PER_TASK - 1 ))  \
    --outputPrefix="/mnt/${out}"       \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05
