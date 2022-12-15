#!/bin/bash

set -eu

dx build 02_saige_variant_test --overwrite

saige_data_dir="/saige_pipeline/data"

sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

output_dir="${saige_data_dir}/02_saige_variant_test/"
dx mkdir --parents "${output_dir}"

# Pheno list
# Obtained using (zcat ukb_obesity_phenos_and_covariates.tsv.gz | head -1 | cut -f31- | sed -e 's/\t/\n/g' | grep -v "invnorm" )
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

brava_phenos=(
  "30690-0.0"
  "30780-0.0"
  "30760-0.0"
  "30870-0.0"
  "30710-0.0"
  "alanine_aminotransferase"
  "30650-0.0"
)

# pheno_col="${phenos[2]}"
pheno_col="${brava_phenos[5]}"

# for i in {0,1,3,4}; do

  # pheno_col="${obesity_phenos[$i]}"
  # pheno_col="${brava_phenos[$i]}"

  # for sex in {both,male,female}; do
sex="both"

if [ "${sex}" = "both" ]; then
  fit_null_output_prefix="${pheno_col}"
elif [ "${sex}" = "male" ] || [ "${sex}" = "female" ]; then
  fit_null_output_prefix="${pheno_col}-${sex}"
else
  echo "Invalid sex: ${sex}" && exit 1
fi

echo "${fit_null_output_prefix}"

model_file="/saige_pipeline/data/01_fit_null/${fit_null_output_prefix}.rda"
variance_ratios="/saige_pipeline/data/01_fit_null/${fit_null_output_prefix}.varianceRatio.txt"

for chrom in {14..16}; do

  if [ ${chrom} -eq 23 ]; then
    chrom="X"
  fi

  plink_bfile="wes_450k:/data/05_export_to_plink/ukb_wes_450k.qced.chr${chrom}"

  output_prefix="saige_variant_test.${fit_null_output_prefix}.chr${chrom}" 
  
  results_file="${output_dir}/${output_prefix}.tsv"

  if [ $( dx ls -l ${results_file} 2> /dev/null | wc -l  ) -eq 0 ]; then
    echo "running chr${chrom}"
    dx run 02_saige_variant_test \
      -iplink_bfile="${plink_bfile}" \
      -ichrom="${chrom}" \
      -isparse_grm="${sparse_grm}" \
      -isparse_grm_samples="${sparse_grm_samples}" \
      -imodel_file="${model_file}" \
      -ivariance_ratios="${variance_ratios}" \
      -ioutput_prefix="${output_prefix}" \
      --name="02_saige_variant_test-${fit_null_output_prefix}-c${chrom}" \
      --destination "${output_dir}" \
      --brief \
      --priority="low" \
      -y
  fi

  sleep 2
done

  # done

# done