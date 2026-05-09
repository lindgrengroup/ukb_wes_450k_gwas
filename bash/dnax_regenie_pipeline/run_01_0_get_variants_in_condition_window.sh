#!/bin/bash

# --- 1. Check Input ---
if [[ -z "$1" ]]; then
    echo "Error: Missing gene symbol."
    echo "Usage: ./script_1_get_window_variants.sh [GENE_SYMBOL]"
    exit 1
fi

GENE_SYMBOL=$1
echo "Processing gene: ${GENE_SYMBOL}"

# --- 2. Define Paths and Parameters ---
IMPUTATION_DIR="/mnt/project/Bulk/Imputation/UKB imputation from genotype"
BGEN_TEMPLATE="${IMPUTATION_DIR}/ukb22828_c%d_b0_v3.bgen"
DESTINATION="/nbaya/outliers/data/imputed_v3_condition_variants/condition_window_variants/"

# Radius of conditioning window 
# (distance in base pairs from start and stop positions of gene)
# 1500000 = 1.5Mb
# +/-1.5Mb Finemapping window from FinnGen paper (Kurki et al. 2023 https://www.nature.com/articles/s41586-022-05473-8)
# RADIUS=1500000
# +/-1Mb window from "Polyfun" in UKB (Weissbrod et al. 2021 https://pmc.ncbi.nlm.nih.gov/articles/PMC7710571/)
# RADIUS=0
RADIUS=1000000
# RADIUS=1500000
# RADIUS=2000000
echo "Window radius: ${RADIUS}"

# --- 3. Find Gene Coordinates (USING HELPER) ---
echo "Fetching coordinates for ${GENE_SYMBOL}..."
gene_coords=$(bash ./_get_gene_coords.sh ${GENE_SYMBOL})

# Check if the helper script failed (exit status not 0) or returned empty
if [[ $? -ne 0 ]] || [[ -z "$gene_coords" ]]; then
    echo "Error: Could not get coordinates for ${GENE_SYMBOL}." >&2
    # The helper script already printed the specific error
    exit 1
fi

# Read the output into variables
read -r CHR START_POS STOP_POS <<< "$gene_coords"
echo "Found ${GENE_SYMBOL} at ${CHR}:${START_POS}-${STOP_POS}"

# --- 4. Prepare BGEN Path ---
BGEN_FILE="$(printf "${BGEN_TEMPLATE}" "$CHR")"
BGEN_FILE_QUOTED=$(printf %q "${BGEN_FILE}")

# --- 5. Define the SAK Job Command ---
# This command will run on the remote worker
read -r -d '' SCRIPT_CMD <<EOF
# Exit on any error
set -e

# --- Part A: Define filenames ---
BGEN_VARIANTS_FILE="imputed_v3_bgen_variants.chr${CHR}.txt"
VARIANTS_IN_WINDOW_FILE="imputed_v3_bgen_variants.${GENE_SYMBOL}.radius${RADIUS}bp.raw.txt"
FINAL_VARIANT_LIST="imputed_v3_bgen_variants.${GENE_SYMBOL}.radius${RADIUS}bp.chrposrefalt.txt"

# --- Part B: Get all variants in the chromosome ---
echo "Extracting all variants for chr${CHR}..."
bgenix -g ${BGEN_FILE_QUOTED} -list > \${BGEN_VARIANTS_FILE}
echo "Extraction complete."

# --- Part C: Filter to window ---
# Use variables expanded from the outer script
WINDOW_START=$(( ${START_POS} - ${RADIUS} ))
WINDOW_STOP=$(( ${STOP_POS} + ${RADIUS} ))

echo "Filtering variants to window ${CHR}:\${WINDOW_START}-\${WINDOW_STOP}..."

# Filter to SNVs (length 1) in the window
awk \
  -v WINDOW_START=\${WINDOW_START} \
  -v WINDOW_STOP=\${WINDOW_STOP} \
  '\$4>=WINDOW_START && \$4 <=WINDOW_STOP && length(\$6) == 1 && length(\$7) == 1' \
  \${BGEN_VARIANTS_FILE} > \${VARIANTS_IN_WINDOW_FILE}

echo "Filtering complete. Reformatting IDs..."

# --- Part D: Reformat IDs to chr:pos:ref:alt ---
# 1st column (alternate_id) is already in format chr:pos_ref_alt
cut -f1 \${VARIANTS_IN_WINDOW_FILE} | sed 's/_/:/g' > \${FINAL_VARIANT_LIST}

# --- Part E: Cleanup ---
rm \${BGEN_VARIANTS_FILE}
rm \${VARIANTS_IN_WINDOW_FILE}

echo "Compressing variant lists..."
gzip \${FINAL_VARIANT_LIST}

echo "Final variant list created: \${FINAL_VARIANT_LIST}.gz"

echo "All steps finished."
EOF

# --- 6. Run the SAK Job ---
priority="high"
instance_type="mem1_ssd1_v2_x8"
name="get_window_variants_radius${RADIUS}bp_${GENE_SYMBOL}_chr${CHR}" 

dx run swiss-army-knife \
  -icmd="$SCRIPT_CMD" \
  --name=$name \
  --instance-type "$instance_type" \
  --priority="$priority" \
  --destination="$DESTINATION" \
  --brief \
  -y

echo "Job submitted. Output will be in ${DESTINATION}"
