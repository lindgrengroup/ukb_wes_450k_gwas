#!/bin/bash

set -eu

dx build 02_saige_variant_test --overwrite

saige_data_dir="/saige_pipeline/data"

sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

output_dir="${saige_data_dir}/02_saige_variant_test/"

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

dx mkdir --parents "${output_dir}"

# for i in {1,2}; do

pheno_col="${phenos[1]}"

fit_null_output_prefix="/saige_pipeline/data/01_fit_null/${pheno_col}"

model_file="${fit_null_output_prefix}.rda"
variance_ratios="${fit_null_output_prefix}.varianceRatio.txt"

# for chrom in {1..23}; do
chrom=21
if [ ${chrom} -eq 23 ]; then
  chrom="X"
fi

plink_bfile="wes_450k:/data/05_export_to_plink/ukb_wes_450k.qced.chr${chrom}"

output_file="saige_variant_test.${pheno_col}.chr${chrom}.tsv"

dx run 02_saige_variant_test \
  -iplink_bfile="${plink_bfile}" \
  -ichrom="${chrom}" \
  -isparse_grm="${sparse_grm}" \
  -isparse_grm_samples="${sparse_grm_samples}" \
  -imodel_file="${model_file}" \
  -ivariance_ratios="${variance_ratios}" \
  -ioutput_file="${output_file}" \
  --name="02_saige_variant_test-${pheno_col}-c${chrom}" \
  --destination "${output_dir}" \
  --brief \
  --priority="low" \
  -y

sleep 2
# done
# done