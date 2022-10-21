#!/bin/bash

saige_data_dir="/saige_pipeline/data"

sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
plink_for_vr_bfile="wes_450k:${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.for_vr"

out_dir="/saige_pipeline/data/01_fit_null/"

# Pheno list:
phenos=(
  whr
  body_mass_index_bmi
  body_fat_percentage
  visceral_adipose_tissue_volume_vat
  abdominal_fat_ratio
  bmi_impedance
  total_tissue_fatp
  android_tissue_fatp
  gynoid_tissue_fatp
  tissuefatp_androidgynoidratio
)
#
# Obtained using (zcat ukb_obesity_phenos_and_covariates.tsv.gz | head -1 | cut -f31- | sed -e 's/\t/\n/g' | grep -v "invnorm" )

pheno_file="${saige_data_dir}/phenotypes/ukb_obesity_phenos_and_covariates.tsv.gz"

dx build 01_saige_fit_null --overwrite
dx mkdir --parents "${out_dir}"

for i in {3..9}; do
  pheno_col="${phenos[$i]}"

  echo $pheno_col
  
  # output_prefix="test_pheno_qt"
  output_prefix="${pheno_col}"

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
done