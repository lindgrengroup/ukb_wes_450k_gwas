#!/bin/bash
#
#

set -eu

saige_data_dir="/saige_pipeline/data"

## OPTIONS
# phenotype_group:
# - "obesity"
# - "biomarkers"
phenotype_group="biomarkers"

## OUTPUT
# Output directory
out_dir="/saige_pipeline/data/01_fit_null/${phenotype_group}/"

# Phenotype list
phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.${phenotype_group}.txt" ) )

sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
plink_for_vr_bfile="wes_450k:${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.for_vr" 

dx build 01_saige_fit_null --overwrite
dx mkdir --parents "${out_dir}"


# pheno_col="${obesity_phenos[3]}"
# pheno_col="alanine_aminotransferase"

for i in {1,2}; do
  pheno_col="${phenos[$i]}"

  sex="both"
  # for sex in {male,female}; do
  # sex="female"

  if [ "${sex}" = "both" ]; then
    
    output_prefix="${pheno_col}"

    pheno_file="${saige_data_dir}/phenotypes/ukb_obesity_phenos_and_covariates.tsv.gz"
    # pheno_file="${saige_data_dir}/phenotypes/ukb_brava_phenos_and_covariates-alt.tsv.gz"

  elif [ "${sex}" = "male" ] || [ "${sex}" = "female" ]; then
    
    output_prefix="${pheno_col}-${sex}"
    
    pheno_file="${saige_data_dir}/phenotypes/ukb_obesity_phenos_and_covariates.${sex}.tsv.gz"
    # pheno_file="${saige_data_dir}/phenotypes/ukb_brava_phenos_and_covariates.${sex}.tsv.gz"

  else
    echo "Invalid sex: ${sex}" && exit 1
  fi
  # output_prefix="test_pheno_qt"

  results_files="${out_dir}/${output_prefix}"

  # Check if output file exists
  if [ $( dx ls -l ${results_files}* 2> /dev/null | wc -l  ) -ne 2 ]; then
    echo "Starting ${output_prefix}"
    dx run 01_saige_fit_null \
      -iplink_for_vr_bfile="${plink_for_vr_bfile}" \
      -isparse_grm="${sparse_grm}" \
      -isparse_grm_samples="${sparse_grm_samples}" \
      -ipheno_file="${pheno_file}" \
      -ipheno_col="${pheno_col}" \
      -itrait_type="${trait_type}" \
      -ioutput_prefix="${output_prefix}" \
      --name="01_saige_fit_null-${output_prefix}" \
      --destination "${out_dir}" \
      --brief \
      --priority="low" \
      -y

  else
    echo "${output_prefix} exists! Skipping."
  fi

  sleep 2

  # done
done