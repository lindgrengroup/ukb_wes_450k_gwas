#!/bin/bash
# 
# Use for running multiple phenotypes for a given chromosome. 
# The walltime is slower than running a completely parallelized set of jobs (one job per phenotype per chrom) but the overall cost should be cheaper because the genotype data for each chrom is only being downloaded once, as opposed to being downloaded for each phenotype in separate jobs.
#

set -eu

dx build 02_saige_all_test_multipheno --overwrite
#echo "WARNING: Did not overwrite 02_saige_all_test_multipheno"

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
# - "longitudinal"
# - "obesity_proteomics"
# - "biomarker_mahalanobis"
# - "mahalanobis_v2_continuous"
# - "prs_as_covariate_v1" # Using genomics PLC 'standard PRS' (trained on non-UKB data)
# - "prs_as_covariate_v2" # In-house PGS that we calculated using PRS-CS applied to UKB imputed data GWAS
# - "prs_as_covariate_v3" # Using genomics PLC 'Enhanced PRS' (informed by UKB data)
# - "mahalanobis_v2_1_irnt" # Covariate-residualized phenotypes, then IRNT-ed
# - "mahalanobis_v2_1_irnt_stdresid" # Covariate-residualized phenotypes, then IRNT-ed, then (genomics PLC standard PRS) regressed out
# - "mahalanobis_v2_1_irnt_lower_vs_inlier" # Outlier classication (cases=lower outlier, controls=inliers)
# - "mahalanobis_v2_1_irnt_upper_vs_inlier" # Outlier classication, Mahalanobis alpha=0.001 (cases=upper outlier, controls=inliers)
# - "mahalanobis_v2_1_irnt_alpha0.01_lower_vs_inlier" # Outlier classication (cases=lower outlier, controls=inliers) with Mahalanobis pval<0.01
# - "mahalanobis_v2_1_irnt_alpha0.01_upper_vs_inlier" # Outlier classication (cases=upper outlier, controls=inliers) with Mahalanobis pval<0.01
# - "original_phenos" # Original phenotypes corresponding to genomics PLC Standard PRS 
# - "standardprs_covariateresid" - standard PRS, stratified by disease status, with covariates residualised out
# - "standardprs_covariateresid_unrelated" - standard PRS, stratified by disease status, with covariates residualised out, individuals 3rd degree or closer (relatedness > 2^(-3.5) removed with greedy algorithm
phenotype_group="biomarker_mahalanobis"

# [INPUT] Sparse GRM
sparse_grm="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_${pop}.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparse_grm_samples="${sparse_grm}.sampleIDs.txt"

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

  # [OPTION] use_irnt:
  # Use SAIGE's built-in inverse-rank normal transformation for quantitative traits
  # - "true" DEFAULT
  # - "false"
  use_irnt="true"

  if [[ "${phenotype_group}" == "mahalanobis_v2_1_"* || "${phenotype_group}" == "standardprs_covariateresid"* ]]; then
    use_irnt="false" # Override to ensure IRNT is not used
  elif [[ "${phenotype_group}" == "biomarker_mahalanobis" ]]; then
    use_irnt="false" # Override to ensure IRNT is not used
  fi

  if [ ${use_irnt} == "true" ]; then
    model_suffix=""
  else
    model_suffix="-no_irnt"
  fi

  # # [OUTPUT] Output directory
  output_dir="${saige_data_dir}/02_saige_all_test${model_suffix}${catevr}/${phenotype_group}/${pop}/${sex}"
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

    echo "running saige_all_test for pheno ids ${pheno_idx_start}-${pheno_idx_stop}, sex=${sex}"
    #for chrom in {1..23}; do
    for chrom in {15..15}; do

      if [ ${chrom} -eq 23 ]; then
        chrom="X"
      fi

      # [INPUT] UKB WES 450k QCed PLINK data
      # plink_bfile="/data/05_export_to_plink/ukb_wes_450k.qced.chr${chrom}" # DEPRECATED v2 (old QC)
      plink_bfile="/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr${chrom}"

      bed_size=$( dx ls -l "${plink_bfile}.bed" 2> /dev/null | cut -f5 -d' ' )
      if (( $( echo "${bed_size} > 240" | bc ) )); then 
        instance_type="mem3_ssd1_v2_x16"
        priority="high"
      elif (( $( echo "${bed_size} > 120" | bc ) )); then 
        instance_type="mem3_ssd1_v2_x8"
        priority="high"
      elif (( $( echo "${bed_size} > 50" | bc ) )); then
        #instance_type="mem3_ssd1_v2_x4"
        instance_type="mem2_ssd1_v2_x8" && echo "############ OVERRIDING TEMPORARILY: mem3_ssd1_v2_x4 machines on low priority are fickle, using mem2_ssd1_v2_x8 instead "
        #priority="low"
        priority="high"
      else
        instance_type="mem3_ssd1_v2_x2"
        priority="low"
      fi

      # TEMPORARY
      #priority="high"
      #echo "###### OVERRIDING TEMPORARILY priority=high ########"

      # [INPUT] Gene consequence annotations
      # group_file="/ukbb-meta/data/annotations/ukb_wes_450k.qced.chr${chrom}.worst_csq_by_gene_canonical.saige.txt.gz" # OLD ANNOTATIONS (no SpliceAI)
      # group_file="/ukb_wes_450k_qc/data/annotations/ukb_wes_450k.qced.brava.v1.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # OLD ANNOTATION (SpliceAI damaging missense not merged)
      # group_file="/ukb_wes_450k_qc/data/annotations/ukb_wes_450k.qced.brava.v2.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Does not include "other_missense" variants
      #group_file="/ukb_wes_450k_qc/data/annotations/brava_v2_csq_w_other_missense/ukb_wes_450k.qced.brava.v2.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Uses variants from old QC
      # group_file="/ukb_wes_450k_qc/data/annotations/brava_v5_csq/ukb_wes_450k.qced.brava.v5.1.saige_group.chr${chrom}.worst_csq_by_gene_canonical.txt.gz" # Only included gnomAD pop max MAF > 0.01 variants
      #group_file="/ukb_wes_450k_qc/data/annotations/brava_v5_csq/ukb_wes_450k.qced.brava.v5.2.saige_group.chr${chrom}.merged.worst_csq_by_gene_canonical.txt.gz" # Old default
      group_file="/ukb_wes_450k_qc/data/annotations/brava_v6_csq/ukb_wes_450k.july.qced.brava.v6.chr${chrom}.tsv.gz" # Current default

      echo "running chr${chrom}"
      dx run 02_saige_all_test_multipheno \
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
        -iuse_irnt="${use_irnt}" \
        --name="02_saige_all_test_multi${model_suffix}-${phenotype_group}-${sex}_${pheno_idx_start}_${pheno_idx_stop}-c${chrom}" \
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
