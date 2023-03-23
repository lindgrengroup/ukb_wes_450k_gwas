#!/usr/bin/env bash
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#
# Author: Nik Baya (2023-03-09)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=saige_02_variant_test
#SBATCH --chdir=/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/saige_step2_variant_test.log
#SBATCH --error=logs/saige_step2_variant_test.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task 4
#SBATCH --requeue
#SBATCH --array=1-23

set -eu # Stop job if any command fails

readonly WD="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"
source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ukb_utils/bash/cluster_utils.sh"

# [OPTION] population
# - "eur"
# - "allpop"
readonly pop="eur"

# [OPTION] phenotype_group:
# - "obesity"
# - "qced_biomarkers"
# - "parkinsons
# - "qced_biomarkers_tail_status_quantile"
# - "qced_biomarkers_tail_status_quantile_midfrac0.683"
# - "smoking"
phenotype_group="obesity"

# [INPUT] Pheno list
if [[ "${phenotype_group}" == *"qced_biomarkers"* ]]; then
  # Always use the same list of biomarkers, with no suffixes related to tail status
  phenos=( $( cat "${WD}/data/phenotypes/phenotype_list.qced_biomarkers.txt" ) )
else
  phenos=( $( cat "${WD}/data/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
fi
readonly n_phenos=${#phenos[@]} # total number of phenotypes in list

# [INPUT] Sparse GRM
readonly sparse_grm="data/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.subset_to_impv3.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
readonly sparse_grm_samples="${sparse_grm}.sampleIDs.txt"

# [OPTION] sex 
# Which sex to run SAIGE on
# - "both_sexes"
# - "female"
# - "male"
# for sex in {male,female}; do
sex="both_sexes"

# Check if sex is valid (i.e. one of "both_sexes", "female", "male")
if [[ ! " ( both_sexes female male ) " =~ " ${sex} " ]]; then
  echo "Invalid sex: ${sex}" && exit 1
fi

# [INPUT] SAIGE null model file directory
fit_null_output_prefix="data/01_fit_null-imputed_v3/${phenotype_group}/${pop}/${sex}"

# [OUTPUT] Output directory
out_dir_singularity="data/02_saige_variant_test-imputed_v3/${phenotype_group}/${pop}/${sex}"
out_dir_bmrc="${WD}/${out_dir_singularity}"
mkdir -p "${out_dir_bmrc}"

# NOTE: pheno_idx is 0-indexed
pheno_idx=$1 # $(( $SLURM_ARRAY_TASK_ID - 1 ))

pheno_col="${phenos[$pheno_idx]}"
gwas_id="${pheno_col}-${pop}-${sex}-imputed_v3"

# [INPUT] SAIGE null model files
model_file="${fit_null_output_prefix}/${gwas_id}.rda"
variance_ratios="${fit_null_output_prefix}/${gwas_id}.varianceRatio.txt"

chrom=$SLURM_ARRAY_TASK_ID

if [ ${chrom} -eq 23 ]; then
  chrom="X"
fi

# [INPUT] Imputed data
#bgen_dir="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_outlier_phens/data/regenie/genotypes/ukb_imputed_v3/minmaf_0.001/info0.8"
#bgen="chr${chrom}.bgen"
#bgen_idx="${bgen}.bgi"
#sample_file="chr${chrom}.sample"
bfile="data/00_set_up/ukb_impv3.subset_to_wes_450k_qc_pass_eur.chr${chrom}"

# [OUTPUT] Output file names
output_prefix="saige_variant_test.${gwas_id}.chr${chrom}" 
results_file_singularity_raw="${out_dir_singularity}/${output_prefix}.tsv"
results_file_bmrc_raw="${out_dir_bmrc}/${output_prefix}.tsv"
results_file_bmrc="${results_file_bmrc_raw}.gz"

# Check if output file exists
if [[ ! -s ${results_file_bmrc} ]]; then

  echo "Starting ${output_prefix}"

  start_time=${SECONDS}
  singularity exec \
    --bind ${WD}:/mnt \
    /apps/singularity/saige_1.1.6.3.sif \
    step2_SPAtests.R  \
      --bedFile "/mnt/${bfile}.bed" \
      --bimFile "/mnt/${bfile}.bim" \
      --famFile "/mnt/${bfile}.fam" \
      --AlleleOrder "ref-first" \
      --minMAF 0  \
      --minMAC 20  \
      --GMMATmodelFile "/mnt/${model_file}"  \
      --varianceRatioFile "/mnt/${variance_ratios}"  \
      --sparseGRMFile "/mnt/${sparse_grm}"  \
      --sparseGRMSampleIDFile "/mnt/${sparse_grm_samples}" \
      --LOCO FALSE  \
      --is_Firth_beta TRUE  \
      --pCutoffforFirth 0.1 \
      --is_output_moreDetails TRUE  \
      --is_fastTest TRUE  \
      --SAIGEOutputFile "/mnt/${results_file_singularity_raw}" \
  && gzip ${results_file_bmrc_raw} \
  && print_update "Finished SAIGE step 2 for ${output_prefix}" $(( ${SECONDS} - ${start_time} ))
else
  print_update "${gwas_id} results exist! Skipping."
fi
