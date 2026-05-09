#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="calc_allele_freq.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}

## OPTIONS
anc='EUR'

# Base pheno group (e.g., prs_as_covariate_v2_2_qt or prs_as_covariate_v2_2_bt)
pheno_group="prs_as_covariate_v2_2_qt"

for pheno_idx in {2..7}; do # Adjust index range as needed
  if [[ "$pheno_group" == "original_phenos_qt" ]]; then
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep -v "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
      current_pheno_group="${pheno_group}"
      
  elif [[ "$pheno_group" == "original_phenos_bt" ]]; then
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
      current_pheno_group="${pheno_group}"
      
  elif [[ "${pheno_group}" == "standardprs_covariateresid_v2_2_ctrls" ]]; then
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.standardprs_covariateresid_v2_2.txt | grep "ctrl" | sed "${pheno_idx}q;d" )
      current_pheno_group="${pheno_group}"
      
  elif [[ "${pheno_group}" == "standardprs_covariateresid_v2_2_cases" ]]; then
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.standardprs_covariateresid_v2_2.txt | grep "case" | sed "${pheno_idx}q;d" )
      current_pheno_group="${pheno_group}"
      
  elif [[ "${pheno_group}" == "prs_as_covariate_v2_2_qt" ]]; then
      # Fetch from standard quantitative trait list
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep -v "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
      # Structure is different: append the trait to the group name
      current_pheno_group="${pheno_group}_${pheno_col}"
      
  elif [[ "${pheno_group}" == "prs_as_covariate_v2_2_bt" ]]; then
      # Fetch from standard binary trait list
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
      # Structure is different: append the trait to the group name
      current_pheno_group="${pheno_group}_${pheno_col}"
      
  else
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.${pheno_group}.txt | sed "${pheno_idx}q;d" )
      current_pheno_group="${pheno_group}"
  fi
  
  # Skip if pheno_col is empty (index out of bounds for the list)
  if [[ -z "${pheno_col}" ]]; then
      continue
  fi

  DESTINATION_DIR="/nbaya/regenie/data/allele_freq" # DNAnexus destination directory

  # Note: Now using $current_pheno_group for all paths and commands
  SAMPLE_ID_PATH="${DESTINATION_DIR}/ukb_wes_450k.qced.EUR.${current_pheno_group}.${pheno_col}.iids.txt"
  if ! dx ls "${SAMPLE_ID_PATH}" &> /dev/null; then
    echo "Sample IDs file does not exist: ${SAMPLE_ID_PATH}. Exiting." && exit 0
  fi
  
  for chrom in {1..23}; do
    if [ $chrom -eq 23 ]; then
      chrom="X"
    fi
  
    OUTPUT_FILE="ukb_wes_450k.qced.chr${chrom}.${anc}.${current_pheno_group}.${pheno_col}.frq"
    FULL_DX_PATH="${DESTINATION_DIR}/${OUTPUT_FILE}"
  
    if ! dx ls "${FULL_DX_PATH}" &> /dev/null; then
    
      echo "Calculating allele frequencies for $anc ${current_pheno_group} $pheno_col for chr$chrom"
      
      dx run swiss-army-knife \
        -iin="${script_dnax}" \
        -icmd="bash ${script} ${anc} ${current_pheno_group} ${pheno_col} ${chrom}" \
        --name="calc_allele_freq_${anc}_${current_pheno_group}_${pheno_col}_chr${chrom}" \
        --instance-type "mem2_ssd1_v2_x8" \
        --priority="high" \
        --destination="${DESTINATION_DIR}" \
        --brief \
        -y
    else
      echo "File ${OUTPUT_FILE} already exists. Skipping chromosome ${chrom}."
    fi
  done
done