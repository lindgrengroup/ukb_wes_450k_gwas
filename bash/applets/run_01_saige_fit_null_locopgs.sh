#!/bin/bash

set -eu

dx build 01_saige_fit_null --overwrite

saige_data_dir="/saige_pipeline/data"

# [OPTION] population:
# - "eur"
# - "allpop"
pop="eur"

# [OPTION] phenotype_group:
# - "locoprs_as_covariate_v2" # Using LOCO PRS-CS on our own GWAS of UKB imputed data as a covariate
phenotype_group="locoprs_as_covariate_v2"

# [OPTION] trait_type:
# - "binary"
# - "quantitative"
trait_type="quantitative"


# [OPTION] use_irnt:
# Use SAIGE's built-in inverse-rank normal transformation for quantitative traits
# - "true" DEFAULT
# - "false" 
use_irnt="true" ## true"


# [INPUT] Sparse GRM
sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${sparse_grm}.sampleIDs.txt"

# [INPUT] PLINK files for variance ratio estimation
plink_for_vr_bfile="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.for_vr" 

# [INPUT] Pheno list
if [[ "${phenotype_group}" == *"qced_biomarkers"* ]]; then
  # Always use the same list of biomarkers, with no suffixes related to tail status
  phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.qced_biomarkers.txt" ) )
else
  phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
fi
n_phenos=${#phenos[@]} # total number of phenotypes in list


# [OPTION] Sex 
# (which sex to run SAIGE on)
# - "both_sexes"
# - "female"
# - "male"
#for sex in {"male","female"}; do
#for sex in {"both_sexes","female","male"}; do
sex="both_sexes"
#sex="female"
  
  
  if [[ ! " ( both_sexes female male ) " =~ " ${sex} " ]]; then
    # Error if sex is invalid (i.e. not one of "both_sexes", "female", "male")
    echo "Invalid sex: ${sex}" && exit 1
  elif [[ "${sex}" == "both_sexes" ]]; then
    # Include sex-related covariates (Sex-combined GWAS)
    sex_covar_col_list=",is_female,is_female_age,is_female_age2"
    sex_qcovar_col_list=",is_female"
  else
    # Exclude sex-combined covariates (Sex-specific GWAS)
    sex_covar_col_list=""
    sex_qcovar_col_list=""
  fi
  
  
  # [OPTION] Covariate columns
  if [[ ${phenotype_group} == "infertility" ]]; then
    covar_col_list="assessment_centre,age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21"
    qcovar_col_list="assessment_centre"
  else
    covar_col_list="age,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15,pc16,pc17,pc18,pc19,pc20,pc21,age2,assessment_centre,sequencing_tranche${sex_covar_col_list}"
    qcovar_col_list="assessment_centre,sequencing_tranche${sex_qcovar_col_list}"
  fi
  
  # [INPUT] Phenotype file
  if [[ ${phenotype_group} == "infertility" ]]; then
    pheno_file="/Laura/phenos/UKBB_WES_Infertility_EUR_ONLY.txt"
  elif [[ ${phenotype_group} == "prs_v2" ]]; then
    pheno_file="${saige_data_dir}/phenotypes/ukb_wes.phenos_and_covariates.prs_as_covariate_v2.tsv.gz"
  else
    # pheno_file="${saige_data_dir}/phenotypes/ukb_${phenotype_group}_phenos_and_covariates.${sex}.tsv.gz"
    # pheno_file="${saige_data_dir}/phenotypes/ukb_wes.${phenotype_group}_phenos_and_covariates.${sex}.tsv.gz"
    pheno_file="${saige_data_dir}/phenotypes/ukb_wes.phenos_and_covariates.${phenotype_group}.tsv.gz"
  fi
  
  # [OUTPUT] Output directory
  if [[ ${use_irnt} == "true" ]]; then
    model_suffix=""
  else
    model_suffix="-no_irnt"
  fi
  out_dir="/saige_pipeline/data/01_fit_null${model_suffix}/${phenotype_group}/${pop}/${sex}"
  dx mkdir --parents "${out_dir}"
  
  # [OPTION] tail_type:
  # Only relevant to tail status phenotype group, otherwise leave as empty string
  # - "" (empty string)
  # - "_qced_is_{lower,upper,either}_tail_quantile{0.01,0.05,0.1}"
  # - "_qced_is_{lower,upper,either}_tail_quantile{0.01,0.05,0.1}_midfrac0.68"
  # for quantile in {0.01,0.1}; do 
  # for tail in {lower,upper}; do
    # tail_type="_qced_is_${tail}_tail_quantile0.01_midfrac0.683"
  tail_type=""
  
  for loco_chrom in {13..22}; do
  
    pheno_col="whr_adj_bmi" #${phenos[$pheno_idx]}${tail_type}"

    # Include PRS as covariate if phenotype group is "prs_as_covariate"
    if [[ ${phenotype_group} == "prs_as_covariate_v1" ]] || [[ ${phenotype_group} == "prs_as_covariate_v3" ]] ; then
      final_covar_col_list="${covar_col_list},prs_${pheno_col}"
    elif [[ ${phenotype_group} == "prs_as_covariate_v2" ]]; then
      final_covar_col_list="${covar_col_list},pgs_${pheno_col}"
    elif [[ ${phenotype_group} == "locoprs_as_covariate_v2" ]]; then
      final_covar_col_list="${covar_col_list},locopgs_${pheno_col}_c${loco_chrom}"
    else
      final_covar_col_list=${covar_col_list}
    fi

    # Only include IRNT-ed residualized phenotypes (specific to Mahalanobis phenotypes)
    #if [[ $pheno_col != *"raw_resid_zscore"* ]]; then 
    #if [[ ${phenotype_group}=="mahalanobis_v2"* ]] && [[ $pheno_col != *"irnt"* ]]; then 
    #  continue
    #fi
  
    # [OUTPUT] Output file names
    gwas_id="${pheno_col}_loco${loco_chrom}-${pop}-${sex}"

    # TEMPORARY
    # Used for testing isLowMemLOCO
    #gwas_id="${gwas_id}-loco"


    results_files="${out_dir}/${gwas_id}"

  
    # Check if output file exists
    if [ $( dx ls -l ${results_files}* 2> /dev/null | wc -l  ) -ne 2 ]; then
    #if [ true ]; then # TEMPORARY - used to save time
      echo "Starting ${gwas_id}"
      dx run 01_saige_fit_null \
        -iplink_for_vr_bed="${plink_for_vr_bfile}.bed" \
        -iplink_for_vr_bim="${plink_for_vr_bfile}.bim" \
        -iplink_for_vr_fam="${plink_for_vr_bfile}.fam" \
        -isparse_grm="${sparse_grm}" \
        -isparse_grm_samples="${sparse_grm_samples}" \
        -ipheno_file="${pheno_file}" \
        -ipheno_col="${pheno_col}" \
        -itrait_type="${trait_type}" \
        -iuse_irnt="${use_irnt}" \
        -isex="${sex}" \
        -icovar_col_list="${final_covar_col_list}" \
        -iqcovar_col_list="${qcovar_col_list}" \
        -ioutput_prefix="${gwas_id}" \
        --name="01_saige_fit_null${model_suffix}-${gwas_id}" \
        --destination "${out_dir}" \
        --brief \
        --priority="low" \
        -y
  
    else
      echo "${gwas_id} exists! Skipping."
    fi
  
    sleep 0.05
  
    # done
  done
  # done
#done
