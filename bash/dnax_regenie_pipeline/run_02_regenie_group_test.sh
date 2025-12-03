#!/bin/bash

set -u # throws error if variables are undefined

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
#readonly pheno_group="mahalanobis_v2_2_irnt_lower_vs_inlier" # Aligned status (lower vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals
#readonly pheno_group="original_phenos_qt" # Quantitative traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="original_phenos_bt" # Dichotomous traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="standardprs_covariateresid_v2_2_ctrls" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
#readonly pheno_group="standardprs_covariateresid_v2_2_cases" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
# readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_ctrls" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche; Also changed the ICD10 code definitions of CAD and osteoporosis
readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_cases" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche; Also changed the ICD10 code definitions of CAD and osteoporosis
#readonly pheno_group="microalbumin_urine_qced" # Only urine microalbumin (QCed following Sinott-Armstrong et al. procedure)

# Variant annotation category
ANNOT_VERSION="v6" # Default used for most of my thesis, including the obesity project
#ANNOT_VERSION="v7" # BRaVa default as of Aug 2025

# Get variables corresponding to `pheno_group`
source _load_variables_for_pheno_group.sh

# Decide whether to use max MAF or max AAF
# - 'MAF' is supposed to replicate SAIGE's cutoffs
# - 'AAF' is the default used by REGENIE
MAF_OR_AAF="MAF"
if [[ "$MAF_OR_AAF" != "MAF" && "$MAF_OR_AAF" != "AAF" ]]; then
  echo "ERROR: MAF_OR_AAF must be 'MAF' or 'AAF'. Got: '$MAF_OR_AAF'" >&2
  exit 1
fi
#AF_CUTOFF=0.001 # Expressed as a fraction (e.g. 0.01 corresponds to 1%)
AF_CUTOFF="0.001" # Expressed as a fraction (e.g. 0.01 corresponds to 1%)
af_flags="--aaf-bins $AF_CUTOFF --vc-maxAAF $AF_CUTOFF"


# --- Conditional Testing Options ---
# condition_flags="" # This variable will no longer be used
CONDITION_LIST_FILE=""
if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
  # Get Ensembl gene ID
  CONDITION_ENSGID="$( bash _get_gene_ensgid.sh "$CONDITION_GENE_IMPUTED" )"
  
  # --- 3. Find Gene Coordinates (USING HELPER) ---
  echo "Fetching coordinates for ${CONDITION_GENE_IMPUTED}..."
  gene_coords=$(bash ./_get_gene_coords.sh ${CONDITION_GENE_IMPUTED})

  if [[ $? -ne 0 ]] || [[ -z "$gene_coords" ]]; then
      echo "Error: Could not get coordinates for ${CONDITION_GENE_IMPUTED}." >&2
      exit 1
  fi

  read -r CONDITION_CHROM_IMPUTED START_POS STOP_POS <<< "$gene_coords"
  echo "Found ${CONDITION_GENE_IMPUTED} at ${CONDITION_CHROM_IMPUTED}:${START_POS}-${STOP_POS}"

  ## Radius of conditioning window, added upstream and downstream of gene start/stop coordinates
  ## Measured in base pairs
#   RADIUS=0 # Minimal
  # RADIUS=1000000 # Default
#   RADIUS=1500000 # Extended
  RADIUS=2000000 # Extended more
  
  # Directory of files of variants lists to condition on
  CONDITION_LIST_DIR="/mnt/project/nbaya/outliers/data/imputed_v3_condition_variants"
  
  ## 1. Define the FULL MOUNTED BGEN path ONCE
  # CONDITION_FILE_MNT_PATH="/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${CONDITION_CHROM_IMPUTED}_b0_v3.bgen"
  CONDITION_FILE_MNT_PATH="${CONDITION_LIST_DIR}/condition_pgen/imputed_v3_qced_snps_maf0.001_hwe1e-10_info0.8_${CONDITION_GENE_IMPUTED}_radius${RADIUS}bp_chr${CONDITION_CHROM_IMPUTED}"
  
  if [[ "${CONDITION_GENE_IMPUTED}" == "HLA-DRB5" ]]; then
    # For this gene, which has a conditioning window that's fully overlapping MHC, use MHC variants to condition, even though Thompson et al. do not use MHC region variants.
    CONDITION_FILE_MNT_PATH+="_with_mhc"
  fi
  
  # 2. Derive all other paths and filenames from it
#   CONDITION_FILE_SAMPLE_MNT_PATH="${CONDITION_FILE_MNT_PATH%.bgen}.sample"
  CONDITION_FILE_FILENAME="${CONDITION_FILE_MNT_PATH##*/}" # This gets just the filename
#   CONDITION_FILE_SAMPLE_FILENAME="${CONDITION_FILE_SAMPLE_MNT_PATH##*/}" # This gets just the filename
  
  # Path to file containing variants (rsid/CHR:POS_REF_ALT if imputed v3; chrCHR:POS:REF:ALT if WES 450k) to condition on (uses --condition-list)
  # CONDITION_LIST_FILE="${CONDITION_LIST_DIR}/imputed_v3_bgen_variants.HMGCR.radius1500000bp.alternate_ids.txt"
  # CONDITION_LIST_FILE="${CONDITION_LIST_DIR}/imputed_v3_bgen_variants.HMGCR.radius1500000bp.rsid.head5000.txt"
  CONDITION_LIST_FILE="tmp-${CONDITION_FILE_MNT_PATH##*/}.txt" # This gets just the filename for the condition pgen prefix and appends ".txt"
  
  
#   echo "WARNING: CONDITION_ENSGID and CONDITION_CHROM_IMPUTED are hard coded for HMGCR"
  
  

fi
# -----------------------------------

# Build flags as an array. This is the most robust way to handle
# arguments that might contain spaces.
flags_array=()

# Add static flags. Assumes flags with values are passed as single items
# (e.g., trait_flag="--trait-file /path/to/trait")
[[ -n "${trait_flag}" ]] && flags_array+=(${trait_flag})
[[ -n "${covar_flag}" ]] && flags_array+=(${covar_flag})
[[ -n "${af_flags}" ]] && flags_array+=(${af_flags})


# If CONDITION_LIST_FILE is not an empty string, build the condition_flags variable
if [[ -n "${CONDITION_LIST_FILE}" ]]; then
  # Add each argument as a separate, quoted element to the array
  flags_array+=(--condition-list="${CONDITION_LIST_FILE}")
  if [[ -n "${CONDITION_FILE_MNT_PATH}" ]]; then
    flags_array+=(
      --condition-file="pgen,${CONDITION_FILE_MNT_PATH}"
      --extract-setlist="${CONDITION_ENSGID}"
      --max-condition-vars=20000
    )
    # --condition-file-sample="${CONDITION_FILE_SAMPLE_FILENAME}"
  fi
fi

PRED="regenie_step1_${anc}_${pheno_group}_pred.list"
STEP1_DIR="nbaya/regenie/data/step1/${anc}"
PRED_PATH="${STEP1_DIR}/${PRED}"
PRED_DNAX="/mnt/project/${STEP1_DIR}/${PRED}"

# pheno_idx: One-indexed (i.e. starts with 1, ends with n)
for pheno_idx in {1..1}; do

  # Get phenotype column name (first column in pred file, in which each row is a different phenotype)
  pheno_col=$( dx cat ${PRED_PATH} | sed "${pheno_idx}q;d" | awk '{ print $1 }')
  
  echo $pheno_col
  pheno_col_flag="--phenoCol=${pheno_col}"
  
  # Create a *copy* of the static flags array for this loop iteration
  local_flags_array=("${flags_array[@]}" "${pheno_col_flag}")
  
  # Convert array to a single string, joined by spaces.
  flags_string="${local_flags_array[*]}"
  
#   flags_string="("${flags_array[@]}" "${pheno_col_flag}")" #$(printf " %q" "${local_flags_array[@]}")
  
  # flags="${trait_flag} ${pheno_col_flag} ${covar_flag} ${af_flags} ${condition_flags}" # This is the old, unsafe method

  LOCO="regenie_step1_${anc}_${pheno_group}_${pheno_idx}.loco"
  LOCO_PATH="${STEP1_DIR}/${LOCO}"

    for chrom in {1..23}; do
      if [ $chrom -eq 23 ]; then
        chrom="X"
      fi

      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        if [ "$chrom" != "$CONDITION_CHROM_IMPUTED" ]; then
          continue
        fi
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
      job_name="regenie_group_test_${ANNOT_VERSION}_${anc}_${pheno_group}_${pheno_col}_c${chrom}"
      
      
      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        out="${out}.conditiongene${CONDITION_GENE_IMPUTED}.radius${RADIUS}bp"
        job_name="${job_name}_radius${RADIUS}bp_cond${CONDITION_GENE_IMPUTED}"
        destination="nbaya/regenie/data/step2/group_tests_conditioned/${pheno_group}/${anc}"
        instance_type="mem3_ssd1_v2_x16"
      else
        destination="nbaya/regenie/data/step2/group_tests/${pheno_group}/${anc}"
      fi
      
      dx mkdir -p ${destination}
      
      # ---  ---
      cmd_prefix=""
      if [[ -n "${CONDITION_LIST_FILE}" ]]; then
        # Create condition list file that includes all variant IDs from the condition pgen
        cmd_prefix="< ${CONDITION_FILE_MNT_PATH}.pvar grep -v '^#CHROM' | cut -f3 > ${CONDITION_LIST_FILE}; " 
      fi
    
      dx run swiss-army-knife \
      	-iin="${script_dnax}" \
        -iin="${bfile_path}.bed" \
        -iin="${bfile_path}.bim" \
        -iin="${bfile_path}.fam" \
        -iin=${LOCO_PATH} \
      	-icmd="${cmd_prefix} bash ${script} ${bfile}.bed ${anc} ${pheno_group} ${PRED_DNAX} ${pheno_col} ${PHENOFILE_DNAX} ${COVARFILE_DNAX} \"${flags_string}\" ${MAF_OR_AAF} ${ANNOT_VERSION} ${chrom} ${out}" \
      	--name="$job_name" \
      	--instance-type "$instance_type" \
      	--priority="$priority" \
      	--destination="$destination" \
      	--brief \
      	-y
    
    done
done
