#!/bin/bash

set -eu

dx build 02_saige_variant_test --overwrite

readonly saige_data_dir="/saige_pipeline/data"

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
# - "infertility"
phenotype_group="qced_biomarkers_tail_status_quantile"

# [OPTION] catevr
# ""        : --isCateVarianceRatio not used for SAIGE step 1
# "-catevr" : --isCateVarianceRatio=TRUE for SAIGE step 1
catevr=""

# [INPUT] Pheno list
if [[ "${phenotype_group}" == *"qced_biomarkers"* ]]; then
  # Always use the same list of biomarkers, with no suffixes related to tail status
  phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.qced_biomarkers.txt" ) )
else
  phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
fi
readonly n_phenos=${#phenos[@]} # total number of phenotypes in list

# [INPUT] Sparse GRM
readonly sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
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
fit_null_output_prefix="${saige_data_dir}/01_fit_null${catevr}/${phenotype_group}/${pop}/${sex}"

# [OUTPUT] Output directory
output_dir="${saige_data_dir}/02_saige_variant_test${catevr}/${phenotype_group}/${pop}/${sex}/"
dx mkdir --parents "${output_dir}"

# [OPTION] tail_type:
# Only relevant to tail status phenotype group, otherwise leave as empty string
# - "" (empty string)
# - "_qced_is_{lower,upper,either}_tail_quantile{0.01,0.05,0.1}"
# for quantile in {0.01,0.1}; do 

for tail in {lower,upper}; do
# tail="lower"
  tail_type="_qced_is_${tail}_tail_quantile0.01"
    #   tail_type=""
  # tail_type=""

  for pheno_idx in {0..32}; do

    pheno_col="${phenos[$pheno_idx]}${tail_type}"  
    gwas_id="${pheno_col}-${pop}-${sex}"

    # [INPUT] SAIGE null model files
    model_file="${fit_null_output_prefix}/${gwas_id}.rda"
    variance_ratios="${fit_null_output_prefix}/${gwas_id}.varianceRatio.txt"
    
    echo "${gwas_id}"

    for chrom in {21..21}; do

      if [ ${chrom} -eq 23 ]; then
        chrom="X"
      fi

      # [INPUT] UKB WES 450k QCed PLINK data
      plink_bfile="/data/05_export_to_plink/ukb_wes_450k.qced.chr${chrom}"

      # [OUTPUT] Output file names
      output_prefix="saige_variant_test${catevr}.${gwas_id}.chr${chrom}" 
      results_file="${output_dir}/${output_prefix}.tsv.gz"

      # Check if output file exists
      if [ $( dx ls -l ${results_file} 2> /dev/null | wc -l  ) -eq 0 ]; then

        echo "Starting chr${chrom}"

        dx run 02_saige_variant_test \
          -ibed="${plink_bfile}.bed" \
          -ibim="${plink_bfile}.bim" \
          -ifam="${plink_bfile}.fam" \
          -isparse_grm="${sparse_grm}" \
          -isparse_grm_samples="${sparse_grm_samples}" \
          -imodel_file="${model_file}" \
          -ivariance_ratios="${variance_ratios}" \
          -ioutput_prefix="${output_prefix}" \
          --name="02_saige_variant_test-${gwas_id}-c${chrom}" \
          --destination "${output_dir}" \
          --brief \
          --priority="low" \
          -y
        
        sleep 0.5 

      fi
       
    done
    
  done

done
  # done
