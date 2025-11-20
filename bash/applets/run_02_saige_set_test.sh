#!/bin/bash

echo "WARNING: use_irnt=false not yet implemented" && exit 1

set -eu

dx build 02_saige_set_test --overwrite

saige_data_dir="/saige_pipeline/data"

# [OPTION] population
# - "eur"
# - "allpop"
readonly pop="eur"

# [OPTION] phenotype_group:
# - "obesity"
# - "qced_biomarkers"
# - "parkinsons
# - "qced_biomarkers_tail_status_quantile"
# - "qced_biomarkers_tail_status_quantile_midfrac0.683"
# - "hormone"
readonly phenotype_group="obesity"

# [OPTION] catevr
# ""        : --isCateVarianceRatio not used for SAIGE step 1
# "-catevr" : --isCateVarianceRatio=TRUE for SAIGE step 1
catevr=""

# [INPUT] Sparse GRM
sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${sparse_grm}.sampleIDs.txt"

# [INPUT] Pheno list
if [[ "${phenotype_group}" == *"qced_biomarkers"* ]]; then
  # Always use the same list of biomarkers, with no suffixes related to tail status
  phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.qced_biomarkers.txt" ) )
else
  phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
fi
n_phenos=${#phenos[@]} # total number of phenotypes in list

# [OPTION] sex 
# Which sex to run SAIGE on
# - "both_sexes"
# - "female"
# - "male"
# for sex in {both_sexes,male,female}; do
# for sex in {male,female}; do
sex="both_sexes"

# Check if sex is valid (i.e. one of "both_sexes", "female", "male")
if [[ ! " ( both_sexes female male ) " =~ " ${sex} " ]]; then
  echo "Invalid sex: ${sex}" && exit 1
fi

# [INPUT] SAIGE null model file directory
fit_null_output_prefix="${saige_data_dir}/01_fit_null${catevr}/${phenotype_group}/${pop}/${sex}"

# [OUTPUT] Output directory
output_dir="${saige_data_dir}/02_saige_set_test${catevr}/${phenotype_group}/${pop}/${sex}"
dx mkdir --parents "${output_dir}"

# [OPTION] tail_type:
# Only relevant to tail status phenotype group, otherwise leave as empty string
# - "" (empty string)
# - "-qced-is_{lower,upper,either}_tail_quantile{0.01,0.05,0.1}"
# for quantile in {0.01,0.1}; do 

#   for tail in {lower,upper,either}; do
    # tail_type="_qced_is_${tail}_tail_quantile${quantile}_midfrac0.683"
tail_type=""

for pheno_idx in `seq 0 ${n_phenos}`; do

  pheno_col="${phenos[$pheno_idx]}${tail_type}" 
  gwas_id="${pheno_col}-${pop}-${sex}"

  # [INPUT] SAIGE null model files
  model_file="${fit_null_output_prefix}/${gwas_id}.rda"
  variance_ratios="${fit_null_output_prefix}/${gwas_id}.varianceRatio.txt"

  echo "${gwas_id}"

  for chrom in {1..23}; do

    if [ ${chrom} -eq 23 ]; then
      chrom="X"
    fi

    # [INPUT] UKB WES 450k QCed PLINK data
    plink_bfile="/data/05_export_to_plink/ukb_wes_450k.qced.chr${chrom}"

    bed_size=$( dx ls -l "${plink_bfile}.bed" 2> /dev/null | cut -f5 -d' ' )
    if (( $( echo "${bed_size} > 120" | bc ) )); then 
      instance_type="mem3_ssd1_v2_x8"
    elif (( $( echo "${bed_size} > 60" | bc ) )); then
      instance_type="mem3_ssd1_v2_x4"
    else
      instance_type="mem3_ssd1_v2_x2"
    fi

    # [INPUT] Gene consequence annotations
    # group_file="/ukbb-meta/data/annotations/ukb_wes_450k.qced.chr${chrom}.worst_csq_by_gene_canonical.saige.txt.gz" # OLD ANNOTATIONS (no SpliceAI)
    # group_file="/ukb_wes_450k_qc/data/annotations/ukb_wes_450k.qced.brava.v1.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # OLD ANNOTATION (SpliceAI damaging missense not merged)
    # group_file="/ukb_wes_450k_qc/data/annotations/ukb_wes_450k.qced.brava.v2.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Does not include "other_missense" variants
    group_file="/ukb_wes_450k_qc/data/annotations/brava_v2_csq_w_other_missense/ukb_wes_450k.qced.brava.v2.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Current Default

    # [OUTPUT] Output files
    output_prefix="saige_set_test${catevr}.${gwas_id}.chr${chrom}"
    results_file="${output_dir}/${output_prefix}.tsv.gz"

    if [ $( dx ls -l ${results_file} 2> /dev/null | wc -l  ) -eq 0 ]; then
      echo "running chr${chrom}"
      dx run 02_saige_set_test \
        -ibed="${plink_bfile}.bed" \
        -ibim="${plink_bfile}.bim" \
        -ifam="${plink_bfile}.fam" \
        -isparse_grm="${sparse_grm}" \
        -isparse_grm_samples="${sparse_grm_samples}" \
        -imodel_file="${model_file}" \
        -ivariance_ratios="${variance_ratios}" \
        -igroup_file="${group_file}" \
        -ioutput_prefix="${output_prefix}" \
        --name="02_saige_set_test-${gwas_id}-c${chrom}" \
        --destination "${output_dir}" \
        --brief \
        --priority="low" \
        --instance-type=${instance_type} \
        -y
    fi

    sleep 1
  done

done

# done
# done
