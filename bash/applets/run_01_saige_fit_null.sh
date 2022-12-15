#!/bin/bash

saige_data_dir="/saige_pipeline/data"

sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
plink_for_vr_bfile="wes_450k:${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.for_vr" 

out_dir="/saige_pipeline/data/01_fit_null/"

# Pheno list
# Columns from ukb_obesity_phenos_and_covariates.tsv.gz
obesity_phenos=(
  body_mass_index_bmi
  whr_adj_bmi
  bmi_impedance
  body_fat_percentage
  visceral_adipose_tissue_volume_vat
  abdominal_fat_ratio
  gynoid_tissue_fatp
  android_tissue_fatp
  total_tissue_fatp
  tissuefatp_androidgynoidratio
)


dx build 01_saige_fit_null --overwrite
dx mkdir --parents "${out_dir}"

pheno_col="${obesity_phenos[3]}"
# pheno_col="alanine_aminotransferase"

# for i in {0..9}; do
#   pheno_col="${phenos[$i]}"

  # for sex in {both,male,female}; do
sex="both"
use_genebass_eur=true

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

if [ ${use_genebass_eur} ]; then
  output_prefix="genebass_eur-${output_prefix}"
fi

echo "${output_prefix}"

dx run 01_saige_fit_null \
  -iplink_for_vr_bfile="${plink_for_vr_bfile}" \
  -isparse_grm="${sparse_grm}" \
  -isparse_grm_samples="${sparse_grm_samples}" \
  -ipheno_file="${pheno_file}" \
  -ipheno_col="${pheno_col}" \
  -ioutput_prefix="${output_prefix}" \
  --name="01_saige_fit_null-${output_prefix}" \
  --destination "${out_dir}" \
  --brief \
  --priority="low" \
  -y

sleep 2

  # done
# done