#!/bin/bash
# 
# Use for running multiple phenotypes for a given chromosome, specifically for genes identified as significant in the obesity project.
#
# The key difference is that the group file this uses is different, it only has annotations for the significant obesity genes.
#
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
# - "qced_biomarkers" - QCed biomarkers
readonly phenotype_group="qced_biomarkers"

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
  if [ ${use_irnt} == "true" ]; then
    model_suffix=""
  else
    model_suffix="-no_irnt"
  fi

  # # [OUTPUT] Output directory
  output_dir="${saige_data_dir}/02_saige_all_test${model_suffix}${catevr}/${phenotype_group}/${pop}/${sex}"
  dx mkdir --parents "${output_dir}"
  
  
  tail_type=""
  
  #for pheno_idx in {0..4}; do
    # 0-indexed
    pheno_idx_start=$1 #${pheno_idx}
    pheno_idx_stop=$2 #${pheno_idx}
    
    echo "running saige_all_test for pheno ids ${pheno_idx_start}-${pheno_idx_stop}, sex=${sex}"
    for chrom in {1..23}; do
    
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
        instance_type="mem2_ssd1_v2_x8"
        priority="low"
      else
        instance_type="mem3_ssd1_v2_x2"
        priority="low"
      fi
    
      # [INPUT] Gene consequence annotations
      group_file="/nbaya/outliers/data/sig_biomarker_mahalanobis_genes_group_file/ukb_wes_450k.july.qced.brava.v6.sig_biomarker_mahalanobis_genes_group_file.Pvalue_Burden_0.001.tsv.gz" # Subset to genes which have Pvalue_Burden < 0.001 in the gene burden test for Mahalanobis distance across 28 QCed biomarkers

    
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
    
      sleep 0.05
    
    done
    #done    
  #done
#done
