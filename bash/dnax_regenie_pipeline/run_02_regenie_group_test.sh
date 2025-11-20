#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="02_regenie_group_test.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}
#dx rm -f ${script_dnax}
#dx upload ${script_local} --path ${script_dnax}

# OPTIONS
anc="EUR" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
#for anc in {"AFR","EAS"}; do

# readonly pheno_group="Height"
#readonly pheno_group="mahalanobis_v2_2_irnt_upper_vs_inlier" # Aligned status (upper vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals
readonly pheno_group="mahalanobis_v2_2_irnt_lower_vs_inlier" # Aligned status (lower vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals
#readonly pheno_group="original_phenos_qt" # Quantitative traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="original_phenos_bt" # Dichotomous traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="standardprs_covariateresid_v2_2_ctrls" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
#readonly pheno_group="standardprs_covariateresid_v2_2_cases" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_ctrls" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche
#readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_cases" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche
#readonly pheno_group="microalbumin_urine_qced" # Only urine microalbumin (QCed following Sinott-Armstrong et al. procedure)

# Variant annotation category
ANNOT_VERSION="v6" # Default used for most of my thesis, including the obesity project
#ANNOT_VERSION="v7" # BRaVa default as of Aug 2025


# Get variables corresponding to pheno_group
export PHENO_GROUP=$pheno_group
source _load_variables_for_pheno_group.sh # This script requires a global variable PHENO_GROUP

# Decide whether to use max MAF or max AAF
MAF_OR_AAF="MAF"
if [[ "$MAF_OR_AAF" != "MAF" && "$MAF_OR_AAF" != "AAF" ]]; then
  echo "ERROR: MAF_OR_AAF must be 'MAF' or 'AAF'. Got: '$MAF_OR_AAF'" >&2
  exit 1
fi
#AF_CUTOFF=0.001 # Expressed as a fraction (e.g. 0.01 corresponds to 1%)
AF_CUTOFF="0.001" # Expressed as a fraction (e.g. 0.01 corresponds to 1%)
af_flags="--aaf-bins $AF_CUTOFF --vc-maxAAF $AF_CUTOFF"

PRED="regenie_step1_${anc}_${pheno_group}_pred.list"
STEP1_DIR="nbaya/regenie/data/step1/${anc}"
PRED_PATH="${STEP1_DIR}/${PRED}"
PRED_DNAX="/mnt/project/${STEP1_DIR}/${PRED}"

#pheno_idx=7 # One-indexed (i.e. starts with 1, ends with n)
for pheno_idx in {2..2}; do

  # Get phenotype column name (first column in pred file, in which each row is a different phenotype)
  pheno_col=$( dx cat ${PRED_PATH} | sed "${pheno_idx}q;d" | awk '{ print $1 }')
  
  echo $pheno_col
  pheno_col_flag="--phenoCol=${pheno_col}"
  
  flags="${trait_flag} ${pheno_col_flag} ${covar_flag} ${af_flags}"
  
    
    LOCO="regenie_step1_${anc}_${pheno_group}_${pheno_idx}.loco"
    LOCO_PATH="${STEP1_DIR}/${LOCO}"
    
    for chrom in {1..1}; do
      if [ $chrom -eq 23 ]; then
        chrom="X"
      fi
    
      # Genotype file
      bfile="ukb_wes_450k.qced.chr${chrom}"
      bfile_path="/Barney/wes/sample_filtered/${bfile}"
    
      echo "Running REGENIE group test with ${bfile} for $anc $pheno_group $pheno_col"
    
      # Allocate machine size based on size of bed file
      bed_size=$( dx ls -l "${bfile_path}.bed" 2> /dev/null | cut -f5 -d' ' )
      if (( $( echo "${bed_size} > 240" | bc ) )); then
        instance_type="mem1_ssd1_v2_x36"
        priority="high"
      elif (( $( echo "${bed_size} > 120" | bc ) )); then
        instance_type="mem1_ssd1_v2_x16"
        priority="high"
      elif (( $( echo "${bed_size} > 50" | bc ) )); then
        instance_type="mem1_ssd1_v2_x8"
        #instance_type="mem2_ssd1_v2_x8" && echo "############ OVERRIDING TEMPORARILY: mem3_ssd1_v2_x4 machines on low priority are fickle, using mem2_ssd1_v2_x8 instead "
        #priority="low"
        priority="high"
      else
        instance_type="mem1_ssd1_v2_x4"
        priority="high"
      fi

      # TEMPORARY OVERRIDE
      #priority="low"
      #if [ $chrom -eq 9 ]; then
      #  priority="high"
      #fi
      priority="high"

      
      out="regenie_group_test.${anc}.${pheno_group}.${pheno_col}.${ANNOT_VERSION}annot.max${MAF_OR_AAF}${AF_CUTOFF}.chr${chrom}"
    
      destination="nbaya/regenie/data/step2/group_tests/${pheno_group}/${anc}"
      dx mkdir -p ${destination}
    
      dx run swiss-army-knife \
      	-iin="${script_dnax}" \
        -iin="${bfile_path}.bed" \
        -iin="${bfile_path}.bim" \
        -iin="${bfile_path}.fam" \
        -iin=${LOCO_PATH} \
      	-icmd="bash ${script} ${bfile}.bed ${anc} ${pheno_group} ${PRED_DNAX} ${pheno_col} ${PHENOFILE_DNAX} ${COVARFILE_DNAX} \"${flags}\" ${MAF_OR_AAF} ${ANNOT_VERSION} ${chrom} ${out}" \
      	--name="regenie_group_test_${ANNOT_VERSION}_${anc}_${pheno_group}_${pheno_col}_c${chrom}" \
      	--instance-type "$instance_type" \
      	--priority="$priority" \
      	--destination="$destination" \
      	--brief \
      	-y
    
    done
done
