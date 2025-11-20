#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="calc_allele_freq_case_ctrl.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}

## OPTIONS
#pheno_group="mahalanobis_v2_2_irnt_upper_vs_inlier"
pheno_group="mahalanobis_v2_2_irnt_lower_vs_inlier"
#pheno_group="original_phenos_qt"

for pheno_idx in {1..7}; do # pheno_idx is 1-indexed (starts at 1, not 0)
	if [[ "$pheno_group" == "original_phenos_qt" ]]; then
      pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep -v "t2d\|cad\|osteoporosis" | sed "${pheno_idx}q;d" )
  else
    pheno_col=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.${pheno_group}.txt | sed "${pheno_idx}q;d" )
  fi

  #pheno_col="Height"
  #pheno_col="microalbumin_urine_qced"
  #pheno_col="t2d_ctrl_prs_resid"
  #pheno_col="aam_bothsexes_raw_resid_irnt_upper_vs_inlier"
  
  DESTINATION_DIR="/nbaya/regenie/data/allele_freq" # DNAnexus destination directory
  
  for chrom in {1..22}; do
    if [ $chrom -eq 23 ]; then
      chrom="X"
    fi
    #readonly BFILE="ukb22418_b0_v2.autosomes" # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")
  
    out="allele_freq_case_ctrl.${pheno_group}.${pheno_col}.chr${chrom}"
    FULL_DX_PATH="${DESTINATION_DIR}/${out}"
  
    if ! dx ls "${FULL_DX_PATH}" &> /dev/null; then
    
      echo "Calculating allele frequencies for $pheno_group $pheno_col for chr$chrom"
      
      dx run swiss-army-knife \
      	-iin="${script_dnax}" \
      	-icmd="pip install pandas; bash ${script} ${pheno_group} ${pheno_col} ${chrom} ${out}" \
      	--name="calc_allele_freq_cc_${pheno_group}_${pheno_col}_chr${chrom}" \
      	--instance-type "mem2_ssd1_v2_x8" \
      	--priority="high" \
      	--destination="${DESTINATION_DIR}" \
      	--brief \
      	-y
    else
      echo "File ${out} already exists. Skipping chromosome ${chrom}."
    fi
  
    
  done
done
