#!/bin/bash

pheno_group=$1
pheno_idx=$2

if [[ "$pheno_group" == "original_phenos_qt" ]]; then
  pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep -v "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
elif [[ "$pheno_group" == "original_phenos_bt" ]]; then
  pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
elif [[ "${pheno_group}" == "standardprs_"*"covariateresid_v2_2_"* ]]; then
  # 1. Get the last part of the string (e.g., "ctrls" or "cases")
  type_suffix="${pheno_group##*_}"

  # 2. Get the base name (e.g., "standardprs_covariateresid_v2_2")
  file_base="${pheno_group%_${type_suffix}}"

  # 3. Get the grep pattern by removing the trailing 's' (e.g., "ctrl" or "case")
  grep_pattern="${type_suffix%s}"

  # 4. Construct the file path
  pheno_list_file="/saige_pipeline/data/phenotypes/phenotype_list.${file_base}.txt"

  pheno_col=$( dx cat "${pheno_list_file}" | grep "${grep_pattern}" | sed "${pheno_idx}q;d" )
#elif [[ "${pheno_group}" == "standardprs_covariateresid_v2_2_cases" ]]; then
#  pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.standardprs_covariateresid_v2_2.txt | grep "case" | sed "${pheno_idx}q;d" )
else
  pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.${pheno_group}.txt | sed "${pheno_idx}q;d" )
fi

echo $pheno_col
