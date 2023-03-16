#!/bin/bash
#
#

set -eu

saige_data_dir="/saige_pipeline/data"

# [OPTION] pop:
# Which genetic ancestry population to subset to.
# - "allpop"
# - "eur"
pop="eur"


# [OPTION] phenotype_group:
# - "obesity"
# - "qced_biomarkers"
phenotype_group="obesity"


# [INPUT] Phenotype list
phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
n_phenos=${#phenos[@]} # total number of phenotypes in list
pheno_file_prefix="${saige_data_dir}/phenotypes/ukb_${phenotype_group}_phenos_and_covariates"

# [INPUT] Sparse GRM
sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${sparse_grm}.sampleIDs.txt"

# [INPUT] PLINK files for variance ratio estimation
plink_for_vr_bfile="wes_450k:${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.for_vr" 

dx build 01_saige_fit_null_quant_trait --overwrite

# [OPTION] sex 
# Which sex to run SAIGE on
# - "both_sexes"
# - "female"
# - "male"
sex="both_sexes"
# for sex in {both_sexes,female,male}; do

# Check if sex is valid (i.e. one of "both_sexes", "female", "male")
if [[ ! " ( both_sexes female male ) " =~ " ${sex} " ]]; then
  echo "Invalid sex: ${sex}" && exit 1
fi

# [INPUT] Phenotype file
pheno_file="${pheno_file_prefix}.${sex}.tsv.gz"

# [OUTPUT] Output directory
out_dir="/saige_pipeline/data/01_fit_null/${phenotype_group}/${pop}/${sex}"
dx mkdir --parents "${out_dir}"


for pheno_idx in {1..1}; do
  pheno_col="${phenos[$pheno_idx]}"

  # [OUTPUT] Output file names
  gwas_id="${pheno_col}-${pop}-${sex}"
  results_files="${out_dir}/${gwas_id}"

  # Check if output file exists
  if [ $( dx ls -l ${results_files}* 2> /dev/null | wc -l  ) -ne 2 ]; then
    
    echo "Starting ${gwas_id}"

    dx run 01_saige_fit_null_quant_trait \
      -iplink_for_vr_bfile="${plink_for_vr_bfile}" \
      -isparse_grm="${sparse_grm}" \
      -isparse_grm_samples="${sparse_grm_samples}" \
      -ipheno_file="${pheno_file}" \
      -ipheno_col="${pheno_col}" \
      -ioutput_prefix="${gwas_id}" \
      --name="01_saige_fit_null_quant_trait-${gwas_id}" \
      --destination "${out_dir}" \
      --brief \
      --priority="low" \
      -y

  else
    echo "${gwas_id} results exist! Skipping."
  fi

  sleep 2

done
# done
# done