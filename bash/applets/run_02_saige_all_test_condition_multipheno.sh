#!/bin/bash

set -eu

dx build 02_saige_all_test_condition_multipheno --overwrite
#echo "WARNING: Did not overwrite 02_saige_all_test_multipheno" # Uncomment this line if not overwriting in the above line

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

# [INPUT] Sparse GRM
sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${sparse_grm}.sampleIDs.txt"


# [INPUT] Sentinel variants
# NOTE: List of all sentinel variants from fine-mapping of common-variant GWAS
all_sentinel_variants="${saige_data_dir}/sentinel_variants/sentinel.log10bf_gt_2.obesity.tsv.gz"

# [OPTION] sex 
# Which sex to run SAIGE on
# - "both_sexes"
# - "female"
# - "male"
#for sex in {both_sexes,male,female}; do
#for sex in {female,male}; do
#sex="female"
#sex="male"
sex="both_sexes"
#for sex in {"both_sexes","female","male"}; do
  
  # Check if sex is valid (i.e. one of "both_sexes", "female", "male")
  if [[ ! " ( both_sexes female male ) " =~ " ${sex} " ]]; then
    echo "Invalid sex: ${sex}" && exit 1
  fi
  
  # [OPTION] catevr
  # ""        : --isCateVarianceRatio not used for SAIGE step 1
  # "-catevr" : --isCateVarianceRatio=TRUE for SAIGE step 1
  catevr=""
      
  # # [OUTPUT] Output directory
  output_dir="${saige_data_dir}/02_saige_all_test_condition${catevr}/${phenotype_group}/${pop}/${sex}"
  dx mkdir --parents "${output_dir}"
  
  # [OPTION] tail_type:
  # Only relevant to tail status phenotype group, otherwise leave as empty string
  # - "" (empty string)
  # - "-qced-is_{lower,upper,either}_tail_quantile{0.01,0.05,0.1}"
  # for quantile in {0.01,0.1}; do 
  
  #   for tail in {lower,upper,either}; do
      # tail_type="_qced_is_${tail}_tail_quantile${quantile}_midfrac0.683"
  tail_type=""
  
  for pheno_idx in {1..1}; do
    # 0-indexed
    pheno_idx_start=${pheno_idx}
    pheno_idx_stop=${pheno_idx}
    
    echo "running saige_all_test_cond for pheno ids ${pheno_idx_start}-${pheno_idx_stop}, sex=${sex}"
    for chrom in {16..16}; do
    
      if [ ${chrom} -eq 23 ]; then
        chrom="X"
      fi
    
      # [INPUT] UKB WES 450k QCed PLINK data
      # plink_bfile="/data/05_export_to_plink/ukb_wes_450k.qced.chr${chrom}" # DEPRECATED v2 (old QC)
      # plink_bfile="/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr${chrom}"# Default when not using merged sentinel variants from imputed data
      plink_bfile="/saige_pipeline/data/sentinel_variants/wes_merged_w_sentinels/ukb_wes_450k.qced.merged_w_imputed_sentinels.chr${chrom}"
    
      bed_size=$( dx ls -l "${plink_bfile}.bed" 2> /dev/null | cut -f5 -d' ' )
      if (( $( echo "${bed_size} > 240" | bc ) )); then 
        instance_type="mem3_ssd1_v2_x32"
        priority="high"
      elif (( $( echo "${bed_size} > 120" | bc ) )); then 
        instance_type="mem3_ssd1_v2_x16"
        priority="high"
      elif (( $( echo "${bed_size} > 50" | bc ) )); then
        instance_type="mem3_ssd1_v2_x8"
        priority="low"
      else
        instance_type="mem3_ssd1_v2_x2"
        priority="low"
      fi
    
      # [INPUT] Gene consequence annotations
      # group_file="/ukbb-meta/data/annotations/ukb_wes_450k.qced.chr${chrom}.worst_csq_by_gene_canonical.saige.txt.gz" # OLD ANNOTATIONS (no SpliceAI)
      # group_file="/ukb_wes_450k_qc/data/annotations/ukb_wes_450k.qced.brava.v1.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # OLD ANNOTATION (SpliceAI damaging missense not merged)
      # group_file="/ukb_wes_450k_qc/data/annotations/ukb_wes_450k.qced.brava.v2.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Does not include "other_missense" variants
      #group_file="/ukb_wes_450k_qc/data/annotations/brava_v2_csq_w_other_missense/ukb_wes_450k.qced.brava.v2.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Uses variants from old QC
      # group_file="/ukb_wes_450k_qc/data/annotations/brava_v5_csq/ukb_wes_450k.qced.brava.v5.1.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Only included gnomAD pop max MAF > 0.01 variants
      #group_file="/ukb_wes_450k_qc/data/annotations/brava_v5_csq/ukb_wes_450k.qced.brava.v5.2.saige_group.chr${chrom}.merged.worst_csq_by_gene_canonical.txt.gz" # Current default
      group_file="/ukb_wes_450k_qc/data/annotations/brava_v6_csq/ukb_wes_450k.july.qced.brava.v6.chr${chrom}.tsv.gz"
    
      echo "running chr${chrom}"
      dx run 02_saige_all_test_condition_multipheno \
        -iphenotype_group="${phenotype_group}" \
        -ipop="${pop}" \
        -isex="${sex}" \
        -icatevr="${catevr}" \
        -ichrom="${chrom}" \
        -ipheno_idx_start=${pheno_idx_start} \
        -ipheno_idx_stop=${pheno_idx_stop} \
        -ibed="${plink_bfile}.bed" \
        -ibim="${plink_bfile}.bim" \
        -ifam="${plink_bfile}.fam" \
        -isparse_grm="${sparse_grm}" \
        -isparse_grm_samples="${sparse_grm_samples}" \
        -igroup_file="${group_file}" \
        -iall_sentinel_variants="${all_sentinel_variants}" \
        --name="02_saige_all_test_cond_multi-${phenotype_group}-${sex}_${pheno_idx_start}_${pheno_idx_stop}-c${chrom}" \
        --destination="${output_dir}" \
        --brief \
        --priority="${priority}" \
        --instance-type="${instance_type}" \
        -y 
    
      sleep 0.2
    done
    #done    
  done
#done
