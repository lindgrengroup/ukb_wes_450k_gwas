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
anc="EUR" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
#for anc in {"AFR","EAS"}; do

# readonly PHENO_GROUP="Height" # For BRaVa
#readonly PHENO_GROUP="mahalanobis_v2_2_irnt_upper_vs_inlier" # Aligned status (upper vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals
#readonly PHENO_GROUP="mahalanobis_v2_2_irnt_lower_vs_inlier" # Aligned status (lower vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals
#readonly PHENO_GROUP="original_phenos_qt" # Quantitative traits (the original phenotypes corresponding to those used in misaligned project)
#readonly PHENO_GROUP="original_phenos_bt" # Dichotomous traits (the original phenotypes corresponding to those used in misaligned project)
#readonly PHENO_GROUP="standardprs_covariateresid_v2_2_ctrls" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
#readonly PHENO_GROUP="standardprs_covariateresid_v2_2_cases" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
# readonly PHENO_GROUP="standardprs_unrelated_covariateresid_v2_2_ctrls" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche; Also changed the ICD10 code definitions of CAD and osteoporosis
#readonly PHENO_GROUP="standardprs_unrelated_covariateresid_v2_2_cases" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche; Also changed the ICD10 code definitions of CAD and osteoporosis
readonly PHENO_GROUP="prs_as_covariate_v2_2_qt" # Subset to unrelateds from v2.2 misaligned pipeline for continuous (quantitative) traits (e.g. height, ldl, bmi)
# readonly pheno_group="prs_as_covariate_v2_2_bt" # Subset to unrelated cases and unrelated controls from v2.2 misaligned pipeline for dichotomous (binary) traits (e.g t2d, cad, osteoporosis)
#readonly PHENO_GROUP="microalbumin_urine_qced" # Only urine microalbumin (QCed following Sinott-Armstrong et al. procedure)


# Variant annotation category
ANNOT_VERSION="v6" # Default used for most of my thesis, including the obesity project
#ANNOT_VERSION="v7" # BRaVa default as of Aug 2025
ANNOT_DIR="/mnt/project/nbaya/regenie/data/annotations/${ANNOT_VERSION}"


# Get variables corresponding to `PHENO_GROUP`
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
CONDITION_LIST_FILE=""

# Condition on imputed v3 variants
if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
  # Options

  # Radius of conditioning window, added upstream and downstream of gene start/stop coordinates
  # Measured in base pairs
  # RADIUS=0 # Minimal
  RADIUS=1000000 # Default
  # RADIUS=1500000 # Extended
  # RADIUS=2000000 # Extended more

  # Whether to use all variants in conditioning window or a subset (e.g. selected by elastic net)
  # Options:
  # - "all" : All variants in conditioning window
  # - "elasticnet_seed<SEED NUMBER>" : Use variants selected by elastic net regression of PRS against imputed variants in window (e.g. "elasticnet_seed1")
  #CONDITION_VARIANT_SUBSET="all"
  CONDITION_VARIANT_SUBSET="elasticnet_seed1"
  # CONDITION_VARIANT_SUBSET="elasticnet_seed2"
  # CONDITION_VARIANT_SUBSET="not_conditioned" # use to run the grouped genes without conditioning

  # Load other variables for conditioning:
  # - CONDITION_ENSGID
  # - CONDITION_CHROM_IMPUTED
  # - CONDITION_FILE_MNT_PATH
  # - CONDITION_LIST_FILE
  source _load_variables_for_conditioning.sh
elif [[ "${PHENO_GROUP}" == "standardprs_unrelated_covariateresid_v2_2"* ]]; then
  # This phenogroup is the only one where we typically condition, thus the note to alert that we aren't conditioning
  echo "Note: Not conditioning on any variants"
fi

# Build flags as an array. This is the most robust way to handle
# arguments that might contain spaces.
flags_array=()

# Add static flags. Assumes flags with values are passed as single items
# (e.g., trait_flag="--trait-file /path/to/trait")
[[ -n "${trait_flag}" ]] && flags_array+=(${trait_flag})
[[ -n "${covar_flag}" ]] && flags_array+=(${covar_flag})
[[ -n "${af_flags}" ]] && flags_array+=(${af_flags})


# If CONDITION_LIST_FILE is not an empty string and CONDITION_VARIANT_SUBSET is not "not_conditioned", build the condition_flags variable
if [[ -n "${CONDITION_LIST_FILE}" ]] && [[ "${CONDITION_VARIANT_SUBSET}" != "not_conditioned" ]]; then
  # Add each argument as a separate, quoted element to the array
  flags_array+=(--condition-list="${CONDITION_LIST_FILE}")
  if [[ -n "${CONDITION_FILE_MNT_PATH}" ]]; then
    flags_array+=(
      --condition-file="pgen,${CONDITION_FILE_MNT_PATH}"
      --extract-setlist="${CONDITION_ENSGID}"
      --max-condition-vars=50000
    )
  fi
fi

PRED="regenie_step1_${anc}_${PHENO_GROUP}_pred.list"
STEP1_DIR="nbaya/regenie/data/step1/${anc}"
PRED_PATH="${STEP1_DIR}/${PRED}"
PRED_DNAX="/mnt/project/${STEP1_DIR}/${PRED}"

# pheno_idx: One-indexed (i.e. starts with 1, ends with n)
for pheno_idx in {2..2}; do

  # Get phenotype column name (first column in pred file, in which each row is a different phenotype)
  pheno_col=$( dx cat ${PRED_PATH} | sed "${pheno_idx}q;d" | awk '{ print $1 }')
  
  echo $pheno_col
  pheno_col_flag="--phenoCol=${pheno_col}"
  
  # Create a *copy* of the static flags array for this loop iteration
  local_flags_array=("${flags_array[@]}" "${pheno_col_flag}")
  
  # Convert array to a single string, joined by spaces.
  flags_string="${local_flags_array[*]}"
  
  LOCO="regenie_step1_${anc}_${PHENO_GROUP}_${pheno_idx}.loco"
  LOCO_PATH="${STEP1_DIR}/${LOCO}"

#     for chrom in {1..22} "X" "GROUPED"; do
    for chrom in {21..21}; do 
      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        # If conditioning, only run the chrom where gene is located
        if [ "$chrom" != "$CONDITION_CHROM_IMPUTED" ]; then
          continue
        fi
      fi

      # Genotype file
      if [[ "${CONDITION_GENE_IMPUTED}" == *"_genes" ]]; then
        # Grouped genes, possibly across chroms
        bfile="${CONDITION_GENE_IMPUTED}_merged"
        bfile_path="/nbaya/outliers/data/ukb_wes_grouped_genes/${bfile}"

        # Define REGENIE annotation and setlist file templates
        ANNOT_TEMPLATE="${ANNOT_DIR}/regenie_FILETYPE.${ANNOT_VERSION}.${CONDITION_GENE_IMPUTED}.txt"
      else
        # All genes in chrom, or
        # Single gene in chrom
        bfile="ukb_wes_450k.qced.chr${chrom}"
        bfile_path="/Barney/wes/sample_filtered/${bfile}"

        # Define REGENIE annotation and setlist file templates
        ANNOT_TEMPLATE="${ANNOT_DIR}/regenie_FILETYPE.${ANNOT_VERSION}.chr${chrom}.txt"
      fi

      echo "Running REGENIE group test with ${bfile} for $anc $PHENO_GROUP $pheno_col"

      # Get DNAnexus machine instance type based on the size of the bed file
      instance_type=$( bash ./_get_instance_type_for_bed_file.sh "${bfile_path}.bed" )
      [[ $instance_type != "mem"* ]] && exit 1 # Exit if instance_type is invalid
      
      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        # Get more memory for conditioning jobs
        instance_type="$( echo $instance_type | sed 's/^mem1/mem2/g' )"
        #instance_type="mem3_ssd1_v2_x32" # Override
      fi
      
      priority="high" # Set job priority (low/normal/high)

      out="regenie_group_test.${anc}.${PHENO_GROUP}.${pheno_col}.${ANNOT_VERSION}annot.max${MAF_OR_AAF}${AF_CUTOFF}.chr${chrom}"
      job_name="regenie_group_test_${ANNOT_VERSION}_${anc}_${PHENO_GROUP}_${pheno_col}_c${chrom}"
      
      if [[ -n "${CONDITION_GENE_IMPUTED}" ]]; then
        destination="nbaya/regenie/data/step2/group_tests_conditioned/${PHENO_GROUP}/${anc}"
        #instance_type="mem3_ssd1_v2_x16"
        out+=".conditiongene${CONDITION_GENE_IMPUTED}.radius${RADIUS}bp"
        job_name="${job_name}_radius${RADIUS}bp_cond${CONDITION_GENE_IMPUTED}"
        
        # If not conditioning on all variants, add suffixes to output and job name
        if [[ "${CONDITION_VARIANT_SUBSET}" != "all" ]]; then
          out+=".${CONDITION_VARIANT_SUBSET}"
          job_name+="_${CONDITION_VARIANT_SUBSET}"
        fi

      else
        destination="nbaya/regenie/data/step2/group_tests/${PHENO_GROUP}/${anc}"
      fi
      
      dx mkdir -p ${destination}
      
      # ---  ---
      cmd_prefix=""
      if [[ -n "${CONDITION_LIST_FILE}" ]] && [[ "${CONDITION_VARIANT_SUBSET}" == "all" ]]; then
        # Create condition list file that includes all variant IDs from the condition pgen
        # NOTE: Only do this if including all variants in conditioning window
        cmd_prefix="< ${CONDITION_FILE_MNT_PATH}.pvar grep -v '^#CHROM' | cut -f3 > ${CONDITION_LIST_FILE}; " 
      fi

      dx run swiss-army-knife \
      	-iin="${script_dnax}" \
        -iin="${bfile_path}.bed" \
        -iin="${bfile_path}.bim" \
        -iin="${bfile_path}.fam" \
        -iin=${LOCO_PATH} \
      	-icmd="${cmd_prefix} bash ${script} ${bfile}.bed ${anc} ${PHENO_GROUP} ${PRED_DNAX} ${pheno_col} ${PHENOFILE_DNAX} ${COVARFILE_DNAX} \"${flags_string}\" ${MAF_OR_AAF} ${ANNOT_TEMPLATE} ${chrom} ${out}" \
      	--name="$job_name" \
      	--instance-type "$instance_type" \
      	--priority="$priority" \
      	--destination="$destination" \
      	--brief \
      	-y
    
    done
done
