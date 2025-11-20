#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="get_sample_list.py"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/python/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}

## OPTIONS
sex="both_sexes"
anc_list='EUR'
#pheno_group="original_phenos_qt"
#pheno_col="aam"
#sex='female'
#pheno_group="standardprs_covariateresid_v2_2_ctrls"
pheno_group="standardprs_covariateresid_v2_2_cases"
#pheno_col="t2d_ctrl_prs_resid"
#pheno_group="mahalanobis_v2_2_irnt_upper_vs_inlier"
#pheno_group="mahalanobis_v2_2_irnt_lower_vs_inlier"

for pheno_idx in {1..3}; do # pheno_idx is 1-indexed (starts at 1, not 0)
  if [[ "$pheno_group" == "original_phenos_qt" ]]; then
    pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep -v "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
  elif [[ "$pheno_group" == "original_phenos_bt" ]]; then
    pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
  elif [[ "${pheno_group}" == "standardprs_covariateresid_v2_2_ctrls" ]]; then
    pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.standardprs_covariateresid_v2_2.txt | grep "ctrl" | sed "${pheno_idx}q;d" )
  elif [[ "${pheno_group}" == "standardprs_covariateresid_v2_2_cases" ]]; then
    pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.standardprs_covariateresid_v2_2.txt | grep "case" | sed "${pheno_idx}q;d" )
  else
    pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.${pheno_group}.txt | sed "${pheno_idx}q;d" )
  fi

  #pheno_col="aam_bothsexes_raw_resid_irnt_upper_vs_inlier"
  
  DESTINATION_DIR="/nbaya/regenie/data/allele_freq" # DNAnexus destination directory
  
  
  #  OUTPUT_FILE="ukb_wes_450k.qced.${anc}.{pheno_col}.iids.txt"
  #  FULL_DX_PATH="${DESTINATION_DIR}/${OUTPUT_FILE}"
  #
  #  if ! dx ls "${FULL_DX_PATH}" &> /dev/null; then
    
      echo "Getting sample list for ${pheno_group} ${pheno_col} ${anc_list}" 
      
      dx run swiss-army-knife \
      	-iin="${script_dnax}" \
      	-icmd="pip install pandas; python3 ${script} ${pheno_group} ${pheno_col} ${anc_list} ${sex}" \
        --name="get_sample_list_${pheno_group}_${pheno_col}_$( echo ${anc_list} | sed 's/,/-/g' )_${sex}" \
      	--instance-type "mem2_ssd1_v2_x4" \
      	--priority="high" \
      	--destination="${DESTINATION_DIR}" \
      	--brief \
      	-y
  #  else
  #    echo "File ${OUTPUT_FILE} already exists. Skipping chromosome ${chrom}."
  #  fi
  #
  
done
