#!/bin/bash

set -u # throws error if variables are undefined

if [[ -n "$@" ]]; then 
  # Parse options
  # Iterate over all arguments passed to the script
  # Use format:
  # ./run_02_regenie_group_test.sh condition_gene_imputed=HMGCR
  for arg in "$@"; do
    case $arg in
      condition_gene_imputed=*)
        CONDITION_GENE_IMPUTED="${arg#*=}"
        ;;
    esac
  done
else
  CONDITION_GENE_IMPUTED=""
fi

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="02_regenie_group_test.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}
#dx rm -f ${script_dnax}
#dx upload ${script_local} --path ${script_dnax}

# OPTIONS
readonly anc="EUR" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)

# Phenotype index range to run (One-indexed).
# - For 'prs_as_covariate_v2_2_*', this controls which LINES to read from phenotype_list.txt
# - For grouped traits, this controls which COLUMNS to read from the pred.list
START_PHENO_IDX=2
END_PHENO_IDX=7

#readonly pheno_group="Height" # Uses BRaVa default covariates.
#readonly pheno_group="mahalanobis_v2_2_irnt_upper_vs_inlier" # Aligned status (upper vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals
#readonly pheno_group="mahalanobis_v2_2_irnt_lower_vs_inlier" # Aligned status (lower vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals 
#readonly pheno_group="original_phenos_qt" # Quantitative traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="original_phenos_bt" # Dichotomous traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="standard_prs_controls" 
# WARNING:standardprs_covariateresid_v2_2 should be run in cases and controls separately - otherwise, the mean imputation in REGENIE step 1 can be messed up due to cases and controls being mutatually exclusive.
#readonly pheno_group="standardprs_covariateresid_v2_2_ctrls" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
#readonly pheno_group="standardprs_covariateresid_v2_2_cases" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
#readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_ctrls" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche
#readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_cases" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche
readonly pheno_group="prs_as_covariate_v2_2_qt" # Subset to unrelateds from v2.2 misaligned pipeline for continuous (quantitative) traits (e.g. height, ldl, bmi)
# readonly pheno_group="prs_as_covariate_v2_2_bt" # Subset to unrelated cases and unrelated controls from v2.2 misaligned pipeline for dichotomous (binary) traits (e.g t2d, cad, osteoporosis)
#readonly pheno_group="microalbumin_urine_qced"

# Variant annotation category
ANNOT_VERSION="v6" # Default used for most of my thesis, including the obesity project
ANNOT_DIR="/mnt/project/nbaya/regenie/data/annotations/${ANNOT_VERSION}"

# Decide whether to use max MAF or max AAF
MAF_OR_AAF="MAF"
if [[ "$MAF_OR_AAF" != "MAF" && "$MAF_OR_AAF" != "AAF" ]]; then
  echo "ERROR: MAF_OR_AAF must be 'MAF' or 'AAF'. Got: '$MAF_OR_AAF'" >&2
  exit 1
fi
AF_CUTOFF="0.001" 
af_flags="--aaf-bins $AF_CUTOFF --vc-maxAAF $AF_CUTOFF"

# --- Main Step 2 Execution Function ---
run_step2() {
  local _pheno_group=$1
  local _anc=$2
  local _start_idx=$3
  local _end_idx=$4
  
  # Get variables corresponding to `PHENO_GROUP`
  export PHENO_GROUP=$_pheno_group
  source _load_variables_for_pheno_group.sh

  # --- Conditional Testing Options ---
  local CONDITION_LIST_FILE=""

  # Condition on imputed v3 variants
  if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
    local RADIUS=1000000 # Default
    local CONDITION_VARIANT_SUBSET="elasticnet_seed1"

    # Load other variables for conditioning
    source _load_variables_for_conditioning.sh
  elif [[ "${PHENO_GROUP}" == "standardprs_unrelated_covariateresid_v2_2"* ]]; then
    echo "Note: Not conditioning on any variants"
  fi

  # Build flags as an array.
  local flags_array=()
  [[ -n "${trait_flag}" ]] && flags_array+=(${trait_flag})
  [[ -n "${covar_flag}" ]] && flags_array+=(${covar_flag})
  [[ -n "${af_flags}" ]] && flags_array+=(${af_flags})

  if [[ -n "${CONDITION_LIST_FILE}" ]] && [[ "${CONDITION_VARIANT_SUBSET}" != "not_conditioned" ]]; then
    flags_array+=(--condition-list="${CONDITION_LIST_FILE}")
    if [[ -n "${CONDITION_FILE_MNT_PATH}" ]]; then
      flags_array+=(
        --condition-file="pgen,${CONDITION_FILE_MNT_PATH}"
        --extract-setlist="${CONDITION_ENSGID}"
        --max-condition-vars=50000
      )
    fi
  fi

  local PRED="regenie_step1_${_anc}_${PHENO_GROUP}_pred.list"
  local STEP1_DIR="nbaya/regenie/data/step1/${_anc}"
  local PRED_PATH="${STEP1_DIR}/${PRED}"
  local PRED_DNAX="/mnt/project/${STEP1_DIR}/${PRED}"

  # Iterate exactly over the range passed into the function
  for pheno_idx in $(seq ${_start_idx} ${_end_idx}); do

    # Get phenotype column name
    local pheno_col=$( dx cat ${PRED_PATH} | sed "${pheno_idx}q;d" | awk '{ print $1 }')
    
    echo "Processing phenotype column: $pheno_col (Index: $pheno_idx)"
    local pheno_col_flag="--phenoCol=${pheno_col}"
    
    local local_flags_array=("${flags_array[@]}" "${pheno_col_flag}")
    local flags_string="${local_flags_array[*]}"
    
    local LOCO="regenie_step1_${_anc}_${PHENO_GROUP}_${pheno_idx}.loco"
    local LOCO_PATH="${STEP1_DIR}/${LOCO}"

    for chrom in {1..22} "X"; do 
      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        if [ "$chrom" != "$CONDITION_CHROM_IMPUTED" ]; then
          continue
        fi
      fi

      local bfile=""
      local bfile_path=""
      local ANNOT_TEMPLATE=""

      # Genotype file
      if [[ "${CONDITION_GENE_IMPUTED}" == *"_genes" ]]; then
        bfile="${CONDITION_GENE_IMPUTED}_merged"
        bfile_path="/nbaya/outliers/data/ukb_wes_grouped_genes/${bfile}"
        ANNOT_TEMPLATE="${ANNOT_DIR}/regenie_FILETYPE.${ANNOT_VERSION}.${CONDITION_GENE_IMPUTED}.txt"
      else
        bfile="ukb_wes_450k.qced.chr${chrom}"
        bfile_path="/Barney/wes/sample_filtered/${bfile}"
        ANNOT_TEMPLATE="${ANNOT_DIR}/regenie_FILETYPE.${ANNOT_VERSION}.chr${chrom}.txt"
      fi

      echo "Running REGENIE group test with ${bfile} for $_anc $PHENO_GROUP $pheno_col"

      local instance_type=$( bash ./_get_instance_type_for_bed_file.sh "${bfile_path}.bed" )
      [[ $instance_type != "mem"* ]] && exit 1 
      
      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        instance_type="$( echo $instance_type | sed 's/^mem1/mem2/g' )"
      fi
      
      local priority="high"

      local out="regenie_group_test.${_anc}.${PHENO_GROUP}.${pheno_col}.${ANNOT_VERSION}annot.max${MAF_OR_AAF}${AF_CUTOFF}.chr${chrom}"
      local job_name="regenie_group_test_${ANNOT_VERSION}_${_anc}_${PHENO_GROUP}_${pheno_col}_c${chrom}"
      local destination=""
      
      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        destination="nbaya/regenie/data/step2/group_tests_conditioned/${PHENO_GROUP}/${_anc}"
        out+=".conditiongene${CONDITION_GENE_IMPUTED}.radius${RADIUS}bp"
        job_name="${job_name}_radius${RADIUS}bp_cond${CONDITION_GENE_IMPUTED}"
        
        if [[ "${CONDITION_VARIANT_SUBSET}" != "all" ]]; then
          out+=".${CONDITION_VARIANT_SUBSET}"
          job_name+="_${CONDITION_VARIANT_SUBSET}"
        fi

      else
        destination="nbaya/regenie/data/step2/group_tests/${PHENO_GROUP}/${_anc}"
      fi
      
      dx mkdir -p ${destination}
      
      local cmd_prefix=""
      if [[ -n "${CONDITION_LIST_FILE}" ]] && [[ "${CONDITION_VARIANT_SUBSET}" == "all" ]]; then
        cmd_prefix="< ${CONDITION_FILE_MNT_PATH}.pvar grep -v '^#CHROM' | cut -f3 > ${CONDITION_LIST_FILE}; " 
      fi

      dx run swiss-army-knife \
        -iin="${script_dnax}" \
        -iin="${bfile_path}.bed" \
        -iin="${bfile_path}.bim" \
        -iin="${bfile_path}.fam" \
        -iin=${LOCO_PATH} \
        -icmd="${cmd_prefix} bash ${script} ${bfile}.bed ${_anc} ${PHENO_GROUP} ${PRED_DNAX} ${pheno_col} ${PHENOFILE_DNAX} ${COVARFILE_DNAX} \"${flags_string}\" ${MAF_OR_AAF} ${ANNOT_TEMPLATE} ${chrom} ${out}" \
        --name="$job_name" \
        --instance-type "$instance_type" \
        --priority="$priority" \
        --destination="$destination" \
        --brief \
        -y
    
    done
  done
}

# --- Pipeline Execution ---
if [[ "${pheno_group}" == "prs_as_covariate_v2_2_"* ]]; then
  # Slices the text file to only get the lines between START_PHENO_IDX and END_PHENO_IDX
  pheno_list=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.$pheno_group.txt | sed -n "${START_PHENO_IDX},${END_PHENO_IDX}p" )
  
  for pheno in $pheno_list; do
    # For these individual pheno runs, the pred.list only has 1 column, so we pass 1 and 1 to the function
    run_step2 "${pheno_group}_${pheno}" "${anc}" 1 1
  done
else
  # For grouped runs, pass the start and end indexes so it loops through the pred.list columns
  run_step2 "$pheno_group" "$anc" "${START_PHENO_IDX}" "${END_PHENO_IDX}"
fi
