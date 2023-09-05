#!/usr/bin/env bash
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
# 
# To run all phenotypes in the phenotype group "obesity":
# sbatch --array=0-9 run_01_saige_fit_null.sh 
#
# Author: Nik Baya (2023-03-09)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=saige_01_fit_null
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/saige_step1_fit_null.log
#SBATCH --error=logs/saige_step1_fit_null.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --nodes 2
#SBATCH --requeue

set -eu # Stop job if any command fails

readonly WD="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"
source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ukb_utils/bash/cluster_utils.sh"

# [OPTION] pop
# Population to use
#   Options:
#   - 'eur': Genetically European
#   - 'allpop': Individuals from all populations who pass QC (i.e. no population filter)
readonly pop="eur"

# [OPTION] phenotype_group:
# - "obesity"
readonly phenotype_group="obesity"

# [OPTION] trait_type:
# - "binary"
# - "quantitative"
readonly trait_type="quantitative"

# [INPUT] Sparse GRM
# readonly sparse_grm="data/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" # DEPRECATED (not subset to imputed v3 samples)
readonly sparse_grm="data/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.subset_to_impv3.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
readonly sparse_grm_samples="${sparse_grm}.sampleIDs.txt"

# [INPUT] PLINK files for variance ratio estimation
readonly plink_for_vr_bfile="data/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.for_vr" 

# [INPUT] Pheno list
if [[ "${phenotype_group}" == *"qced_biomarkers"* ]]; then
  # Always use the same list of biomarkers, with no suffixes related to tail status
  readonly phenos=( $( cat "${WD}/data/phenotypes/phenotype_list.qced_biomarkers.txt" ) )
else
  readonly phenos=( $( cat "${WD}/data//phenotypes/phenotype_list.${phenotype_group}.txt" ) )
fi
readonly n_phenos=${#phenos[@]} # total number of phenotypes in list


# [OPTION] Sex 
# (which sex to run SAIGE on)
# - "both_sexes"
# - "female"
# - "male"
sex="female"

# Check if sex is valid (i.e. one of "both_sexes", "female", "male")
if [[ ! " ( both_sexes female male ) " =~ " ${sex} " ]]; then
  echo "Invalid sex: ${sex}" && exit 1
fi

# [INPUT] Phenotype file
pheno_file="data/phenotypes/ukb_imputed_v3_${phenotype_group}_phenos_and_covariates.${sex}.tsv"

# [OUTPUT] Output directory
out_dir_singularity="data/01_fit_null-imputed_v3/${phenotype_group}/${pop}/${sex}" # directory for singularity container
out_dir_bmrc="${WD}/${out_dir_singularity}" # path in BMRC
mkdir -p "${out_dir_bmrc}"

# Set up flags
if [[ "${sex}" == "male" ]] || [[ "${sex}" == "female" ]]; then
  # Do not include sex-related covariates if running sex-specific analysis
  sex_covar_col_list=""
  sex_qcovar_col_list=""
else
  # Include sex-related covariates if not running sex-specific analysis
  sex_covar_col_list=",is_female,is_female_age,is_female_age2"
  sex_qcovar_col_list=",is_female"
fi

# Get number of threads
n_threads=$(( $( nproc --all ) -1 ))

# Get inverse-normalize flag if trait_type=="quantitative"
if [ ${trait_type} == "quantitative" ]; then
  trait_flags="--traitType=${trait_type}   --invNormalize=TRUE"
else
  trait_flags="--traitType=${trait_type}"
fi


# NOTE: pheno_idx is 0-indexed
pheno_idx="${SLURM_ARRAY_TASK_ID}"
pheno_col="${phenos[$pheno_idx]}"

# [OUTPUT] Output file names
gwas_id="${pheno_col}-${pop}-${sex}-imputed_v3"
results_files_singularity="${out_dir_singularity}/${gwas_id}" # path for singularity container
results_files_bmrc="${out_dir_bmrc}/${gwas_id}" # path in BMRC

# Check if output file exists
if [[ ! -s ${results_files_bmrc}.rda ]] || [[ ! -s ${results_files_bmrc}.varianceRatio.txt ]]; then
  rm -f ${results_files_bmrc}.*
  start_time=${SECONDS}
  singularity exec \
    --bind ${WD}:/mnt \
    /apps/singularity/saige_1.1.6.3.sif \
    step1_fitNULLGLMM.R  \
    --plinkFile "/mnt/${plink_for_vr_bfile}" \
    --sparseGRMFile "/mnt/${sparse_grm}" \
    --sparseGRMSampleIDFile "/mnt/${sparse_grm_samples}"  \
    --useSparseGRMtoFitNULL TRUE  \
    --phenoFile "/mnt/${pheno_file}" \
    --skipVarianceRatioEstimation FALSE \
    --isCateVarianceRatio FALSE \
    --IsOverwriteVarianceRatioFile TRUE \
    --phenoCol "${pheno_col}" \
    --covarColList="age,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15,pc16,pc17,pc18,pc19,pc20,pc21,age2,assessment_centre,genotyping_batch${sex_covar_col_list}" \
    --qCovarColList="assessment_centre,genotyping_batch${sex_qcovar_col_list}"  \
    --sampleIDColinphenoFile="IID" \
    ${trait_flags} \
    --outputPrefix="/mnt/${results_files_singularity}" \
    --nThreads=${n_threads}
  if [ $? -eq 0 ]; then
    print_update "Finished SAIGE step 1 for ${gwas_id}" $(( ${SECONDS} - ${start_time} ))
  else
    raise_error "SAIGE step 1 failed for ${gwas_id}"
  fi
else
  echo "${gwas_id} exists! Skipping."
fi
