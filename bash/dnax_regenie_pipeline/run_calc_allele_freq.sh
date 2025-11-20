#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="calc_allele_freq.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}

## OPTIONS
#anc="SAS" # Genetic ancestry group (options: "EUR", "AFR", "EAS", "SAS", "AMR")
anc='EUR'
#for anc in {"EAS","SAS","AMR"}; do
#pheno_group="mahalanobis_v2_2_irnt_upper_vs_inlier"
#pheno_group="mahalanobis_v2_2_irnt_lower_vs_inlier"
#pheno_group="original_phenos_qt"
#pheno_group="original_phenos_bt"
#pheno_group="standardprs_covariateresid_v2_2_ctrls"
pheno_group="standardprs_covariateresid_v2_2_cases"

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
  

  #pheno_col="Height"
  #pheno_col="microalbumin_urine_qced"
  #pheno_col="t2d_ctrl_prs_resid"
  #pheno_col="aam_bothsexes_raw_resid_irnt_upper_vs_inlier"
  
  DESTINATION_DIR="/nbaya/regenie/data/allele_freq" # DNAnexus destination directory

  SAMPLE_ID_PATH="${DESTINATION_DIR}/ukb_wes_450k.qced.EUR.${pheno_group}.${pheno_col}.iids.txt"
  if ! dx ls "${SAMPLE_ID_PATH}" &> /dev/null; then
    echo "Sample IDs file does not exist: ${SAMPLE_ID_PATH}. Exiting." && exit 0
  fi
  
  for chrom in {1..23}; do
    if [ $chrom -eq 23 ]; then
      chrom="X"
    fi
    #readonly BFILE="ukb22418_b0_v2.autosomes" # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")
  
    OUTPUT_FILE="ukb_wes_450k.qced.chr${chrom}.${anc}.${pheno_group}.${pheno_col}.frq"
    #OUTPUT_FILE="ukb_wes_450k.qced.chr${chrom}.${anc}.${pheno_group}.${pheno_col}.frq.count"
    FULL_DX_PATH="${DESTINATION_DIR}/${OUTPUT_FILE}"
  
    if ! dx ls "${FULL_DX_PATH}" &> /dev/null; then
    
      echo "Calculating allele frequencies for $anc ${pheno_group} $pheno_col for chr$chrom"
      
      dx run swiss-army-knife \
      	-iin="${script_dnax}" \
      	-icmd="bash ${script} ${anc} ${pheno_group} ${pheno_col} ${chrom}" \
      	--name="calc_allele_freq_${anc}_${pheno_group}_${pheno_col}_chr${chrom}" \
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
