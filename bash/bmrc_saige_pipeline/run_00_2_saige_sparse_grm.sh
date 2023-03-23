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
#SBATCH --nodes=20
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --requeue

set -e # Stop job if any command fails

readonly WD="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"
source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ukb_utils/bash/cluster_utils.sh"
readonly nnodes=$SLURM_JOB_NUM_NODES
readonly nthreads=200 #$(( $nnodes * 5 )) # $(( $( nproc --all ) - 1 ))
readonly mem_avail_mb=$(( $( grep "MemTotal" /proc/meminfo | cut -f 8 -d' ' ) / 10000 )) # memory available (Mb)
echo "memory available (Mb): $mem_avail_mb"

# [OPTION] pop
# Population to use
#   Options:
#   - 'eur': Genetically European
#   - 'allpop': Individuals from all populations who pass QC (i.e. no population filter)
readonly pop="eur"

# [INPUT]
readonly bfile="ukb_array.wes_450k_qc_pass_${pop}.subset_to_impv3.pruned"

# [OUTPUT]
readonly out_dir="${WD}/data/00_set_up"
readonly out="short-${nnodes}nodes-${nthreads}threads-ukb_array.wes_450k_qc_pass_${pop}.subset_to_impv3.pruned"

mkdir -p ${out_dir}

if [[ ! -f "${out_dir}/${out}.mtx" ]]; then
  start_time=${SECONDS}
  singularity exec \
    --bind ${WD}/data/00_set_up:/mnt \
    /apps/singularity/saige_1.1.6.3.sif \
    createSparseGRM.R  \
      --plinkFile="/mnt/${bfile}" \
      --nThreads=${nthreads}  \
      --outputPrefix="/mnt/${out}"       \
      --numRandomMarkerforSparseKin=5000      \
      --relatednessCutoff=0.05

  if [ $? -eq 0 ]; then
    print_update "Finished creating sparse GRM" $(( ${SECONDS} - ${start_time} ))
  fi
fi
