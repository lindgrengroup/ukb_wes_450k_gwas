#!/bin/bash

# Directory of files of variants lists to condition on
CONDITION_LIST_DIR="/mnt/project/nbaya/outliers/data/imputed_v3_condition_variants"

# Condition file prefix
CONDITION_FILE_PREFIX="imputed_v3_qced_snps_maf0.001_hwe1e-10_info0.8_${CONDITION_GENE_IMPUTED}_radius${RADIUS}bp"

# Get Ensembl gene ID (default GRCh version: 38, because UKB WES uses build 38)
# This ENSGID is used to subset the burden testing in UKB WES to just the single gene
# If grouped gene, returns CONDITION_GENE_IMPUTED unchanged
CONDITION_ENSGID="$( bash _get_gene_ensgid.sh "$CONDITION_GENE_IMPUTED" "38" )"

if [[ "${CONDITION_GENE_IMPUTED}" == *"_genes" ]]; then 
  CONDITION_CHROM_IMPUTED="GROUPED"
else
  # Single gene
  # --- 3. Find Gene Coordinates (USING HELPER) ---
  echo "Fetching coordinates for ${CONDITION_GENE_IMPUTED}..."
  gene_coords=$(bash ./_get_gene_coords.sh ${CONDITION_GENE_IMPUTED})

  if [[ $? -ne 0 ]] || [[ -z "$gene_coords" ]]; then
    echo "Error: Could not get coordinates for ${CONDITION_GENE_IMPUTED}." >&2
    exit 1
  fi

  read -r CONDITION_CHROM_IMPUTED START_POS STOP_POS <<< "$gene_coords"
  echo "Found ${CONDITION_GENE_IMPUTED} at ${CONDITION_CHROM_IMPUTED}:${START_POS}-${STOP_POS}"

  # Add chromosome to end of file path
  CONDITION_FILE_PREFIX+="_chr${CONDITION_CHROM_IMPUTED}"
fi

## 1. Define the FULL MOUNTED BGEN path ONCE
# CONDITION_FILE_MNT_PATH="/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${CONDITION_CHROM_IMPUTED}_b0_v3.bgen"
CONDITION_FILE_MNT_PATH="${CONDITION_LIST_DIR}/condition_pgen/${CONDITION_FILE_PREFIX}"

if [[ "${CONDITION_GENE_IMPUTED}" == "HLA-DRB5" ]]; then
  # For this gene, which has a conditioning window that's fully overlapping MHC, use MHC variants to condition, even though Thompson et al. do not use MHC region variants.
  CONDITION_FILE_MNT_PATH+="_with_mhc"
fi

# 2. Derive all other paths and filenames from it
#   CONDITION_FILE_SAMPLE_MNT_PATH="${CONDITION_FILE_MNT_PATH%.bgen}.sample"
CONDITION_FILE_FILENAME="${CONDITION_FILE_MNT_PATH##*/}" # This gets just the filename
#   CONDITION_FILE_SAMPLE_FILENAME="${CONDITION_FILE_SAMPLE_MNT_PATH##*/}" # This gets just the filename

# Path to file containing variants (rsid/CHR:POS_REF_ALT if imputed v3; chrCHR:POS:REF:ALT if WES 450k) to condition on (uses --condition-list)
if [[ "${CONDITION_VARIANT_SUBSET}" == "all" ]]; then
  CONDITION_LIST_FILE="tmp-${CONDITION_FILE_MNT_PATH##*/}.txt" # This gets just the filename for the condition pgen prefix and appends ".txt"
elif [[ "${CONDITION_VARIANT_SUBSET}" == "elasticnet"* ]]; then
  # Extract seed from CONDITION_VARIANT_SUBSET
  ELASTIC_NET_SEED="${CONDITION_VARIANT_SUBSET##*seed}"

  # Use list of variants selected by elastic net regression (PRS ~ imputed variants in conditioning window)
  CONDITION_LIST_FILE="${CONDITION_LIST_DIR}/lasso_elastic_net/${CONDITION_FILE_PREFIX}_selected_variants_seed${ELASTIC_NET_SEED}.txt"
else
  echo "Error: Invalid CONDITION_VARIANT_SUBSET: ${CONDITION_VARIANT_SUBSET}" >&2
  exit
fi
