#!/bin/bash
# NOTE: This requires running run_get_variants_in_condition_window.sh first, to get the variant lists for each gene conditioning window (for a given RADIUS of conditioning, e.g. RADIUS=1000000bp)


# --- 1. Check Input ---
if [[ -z "$1" ]]; then
    echo "Error: Missing gene symbol."
    echo "Usage: ./run_filter_bgen_to_pgen_fast.sh [GENE_SYMBOL]"
    exit 1
fi

GENE_SYMBOL=$1
echo "Processing gene: ${GENE_SYMBOL}"

# --- 2. Define Paths and Parameters ---
# REMOVED trailing slash
DNAX_WD="/nbaya/outliers/data/imputed_v3_condition_variants"
DESTINATION="${DNAX_WD}/condition_pgen/"
IMPUTATION_DIR="/mnt/project/Bulk/Imputation/UKB imputation from genotype"
BGEN_TEMPLATE="${IMPUTATION_DIR}/ukb22828_c%s_b0_v3.bgen"
SAMPLE_TEMPLATE="${IMPUTATION_DIR}/ukb22828_c%s_b0_v3.sample"
MFI_TEMPLATE="${IMPUTATION_DIR}/ukb22828_c%s_b0_v3.mfi.txt"
FAM_KEEP_LIST="/mnt/project/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr1.fam" # Can be any chromosome

priority="high"
instance_type="mem1_ssd1_v2_x16"

MAF_THRESHOLD="0.001"
HWE_THRESHOLD="1e-10"
INFO_THRESHOLD="0.8"
RADIUS=0 # Minimal
# RADIUS=1000000 # Default
# RADIUS=1500000 # Extended
# RADIUS=2000000 # Extended more

MAF_LABEL="maf$(echo ${MAF_THRESHOLD} | sed 's/\\./p/')"
HWE_LABEL="hwe${HWE_THRESHOLD}"
INFO_LABEL="info$(echo ${INFO_THRESHOLD} | sed 's/\\./p/')"

# --- 3. Find Gene Coordinates (USING HELPER) ---
echo "Fetching coordinates for ${GENE_SYMBOL}..."
gene_coords=$(bash ./_get_gene_coords.sh ${GENE_SYMBOL})

if [[ $? -ne 0 ]] || [[ -z "$gene_coords" ]]; then
    echo "Error: Could not get coordinates for ${GENE_SYMBOL}." >&2
    exit 1
fi

read -r CHR START_POS STOP_POS <<< "$gene_coords"
echo "Found ${GENE_SYMBOL} at ${CHR}:${START_POS}-${STOP_POS}"

# --- 4. Build Paths for this Job ---
CURRENT_BGEN_FILE=$(printf "${BGEN_TEMPLATE}" "$CHR")
CURRENT_SAMPLE_FILE=$(printf "${SAMPLE_TEMPLATE}" "$CHR")
CURRENT_MFI_FILE=$(printf "${MFI_TEMPLATE}" "$CHR")

# Mounted gzipped variant list path
VARIANT_LIST_MNT_PATH="/mnt/project${DNAX_WD}/condition_window_variants/imputed_v3_bgen_variants.${GENE_SYMBOL}.radius${RADIUS}bp.chrposrefalt.txt.gz"
OUTPUT_PREFIX="imputed_v3_qced_snps_${MAF_LABEL}_${HWE_LABEL}_${INFO_LABEL}_${GENE_SYMBOL}_radius${RADIUS}bp_chr${CHR}"
JOB_NAME="filter_pgen_muyrápido_radius${RADIUS}bp_${GENE_SYMBOL}_chr${CHR}"

# --- 5. Define the SAK Job Command ---
read -r -d '' SCRIPT_CMD <<EOF

set -e
echo "plink2 version: \$(plink2 --version)"

# Define Inputs
ORIGINAL_FAM_LIST="${FAM_KEEP_LIST}"
MFI_FILE_PATH="${CURRENT_MFI_FILE}"
BGEN_PATH="${CURRENT_BGEN_FILE}"
SAMPLE_PATH="${CURRENT_SAMPLE_FILE}"

# Define Temps
CORRECTED_FAM_PATH="wes_450k_keep_list_corrected.fam"
INFO_FAILED_LIST="info_failed_variants.snplist"
TEMP_SMALL_BGEN="temp_window.bgen"

# QC Vars
MAF_THRESH=${MAF_THRESHOLD}
HWE_THRESH=${HWE_THRESHOLD}
INFO_THRESH=${INFO_THRESHOLD}

# Calculate Window for bgenix (Gene +/- Radius)
# We use 0 as min start to prevent negative numbers
RAW_START=$(( ${START_POS} - ${RADIUS} ))
WINDOW_START=\$(( RAW_START > 0 ? RAW_START : 0 ))
WINDOW_STOP=$(( ${STOP_POS} + ${RADIUS} ))

# --- LOGIC: Pad Chromosome ONLY for bgenix internal query ---
if [[ "$CHR" =~ ^[1-9]$ ]]; then
    CHR_PADDED="0$CHR"
else
    CHR_PADDED="$CHR"
fi

# This variable will hold the range(s) for bgenix -incl-range
BGENIX_RANGE_CMD="\${CHR_PADDED}:\${WINDOW_START}-\${WINDOW_STOP}"
OUTPUT_MHC_FLAG="" # Flag for output prefix modification

echo "--- Step 1: Pre-filtering BGEN with bgenix ---"
echo "Targeting Gene Window: \${CHR_PADDED}:\${WINDOW_START}-\${WINDOW_STOP}"
echo "BGEN File: \${BGEN_PATH}"

# MHC/HLA Exclusion Configuration (GRCh37/hg19)
# Region coordinates from Gourraud et al. 2014 https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0097282
# 6:28,866,528–33,775,446 (GRCh37/hg19 coordinates)
MHC_CHR_PADDED="06"
MHC_START="28866528"
MHC_END="33775446"

if [ "\$CHR_PADDED" == "\$MHC_CHR_PADDED" ]; then
    
    # 1. Calculate the intersection of the Gene Window and the Pre-MHC region
    PRE_MHC_END=\$(( WINDOW_STOP < MHC_START ? WINDOW_STOP : MHC_START - 1 ))
    
    # 2. Calculate the intersection of the Gene Window and the Post-MHC region
    POST_MHC_START=\$(( WINDOW_START > MHC_END ? WINDOW_START : MHC_END + 1 ))

    # Create the range string(s)
    TEMP_RANGE_CMD=""
    
    # Check if the Pre-MHC or Post-MHC range is non-empty
    IS_PARTIAL_OVERLAP=false
    
    # Add Pre-MHC range if it overlaps with the window AND is valid
    if [ "\$WINDOW_START" -le "\$PRE_MHC_END" ] && [ "\$PRE_MHC_END" -ge 0 ]; then
        TEMP_RANGE_CMD="\${CHR_PADDED}:\${WINDOW_START}-\${PRE_MHC_END}"
        IS_PARTIAL_OVERLAP=true
    fi

    # Add Post-MHC range if it overlaps with the window AND is valid
    if [ "\$POST_MHC_START" -le "\$WINDOW_STOP" ]; then
        if [ -n "\$TEMP_RANGE_CMD" ]; then
            TEMP_RANGE_CMD="\${TEMP_RANGE_CMD},\${CHR_PADDED}:\${POST_MHC_START}-\${WINDOW_STOP}"
        else
            TEMP_RANGE_CMD="\${CHR_PADDED}:\${POST_MHC_START}-\${WINDOW_STOP}"
        fi
        IS_PARTIAL_OVERLAP=true
    fi

    if [ "\$IS_PARTIAL_OVERLAP" == "true" ]; then
        # Use the smart, non-contiguous range that excludes MHC
        BGENIX_RANGE_CMD="\$TEMP_RANGE_CMD"
        echo "Excluding MHC region, using ranges: \$BGENIX_RANGE_CMD"
    else
        # Gene window fully overlaps MHC (IS_PARTIAL_OVERLAP is false and it's Chr 6)
        # Revert to the full window and add the MHC flag to the output prefix
        OUTPUT_MHC_FLAG="_with_mhc"
        BGENIX_RANGE_CMD="\${CHR_PADDED}:\${WINDOW_START}-\${WINDOW_STOP}"
        echo "Warning: Gene window entirely overlaps the MHC region. Including MHC variants and setting output flag: \$OUTPUT_MHC_FLAG"
    fi

fi

# Re-define OUTPUT_PREFIX to include the flag
OUTPUT_PREFIX_NEW="${OUTPUT_PREFIX}\${OUTPUT_MHC_FLAG}"
echo "Final Output Prefix: \${OUTPUT_PREFIX_NEW}"


# Check for index path
BGI_PATH="\${BGEN_PATH}.bgi"
if [ ! -f "\$BGI_PATH" ]; then
    echo "Error: Explicit .bgi index file not found at \${BGI_PATH}" >&2
    exit 1
fi

# Run bgenix with PADDED chromosome and the selected inclusive range(s)
bgenix -g "\${BGEN_PATH}" -i "\$BGI_PATH" -incl-range "\$BGENIX_RANGE_CMD" > "\${TEMP_SMALL_BGEN}"

echo "Created temporary small BGEN file."

echo "--- Step 2: Prep QC Files ---"

# Check variant list
if [ ! -f "${VARIANT_LIST_MNT_PATH}" ]; then
    echo "Error: Variant list file not found at ${VARIANT_LIST_MNT_PATH}"
    exit 1
fi
count_keep_variants=\$((\$(gunzip -c ${VARIANT_LIST_MNT_PATH} | wc -l)))
echo "Found \${count_keep_variants} variants in allow-list."

# INFO Score Filter
echo "Filtering MFI by INFO < \${INFO_THRESH}..."
awk -v chr="$CHR" -v min_info=\${INFO_THRESH} '\$8 < min_info {print chr":"\$3":"\$4":"\$5}' "\${MFI_FILE_PATH}" > \${INFO_FAILED_LIST}

# Fix FAM file
echo "Correcting FAM file..."
awk '{print \$2, \$2, \$3, \$4, \$5, \$6}' "\${ORIGINAL_FAM_LIST}" > "\${CORRECTED_FAM_PATH}"


echo "--- Step 3: Running PLINK2 on small BGEN ---"
# Note: We point --bgen to the TEMP_SMALL_BGEN
# Use the modified OUTPUT_PREFIX_NEW

plink2 \
    --bgen "\${TEMP_SMALL_BGEN}" ref-first \
    --sample "\${SAMPLE_PATH}" \
    --keep "\${CORRECTED_FAM_PATH}" \
    --set-all-var-ids "@:#:\\\$r:\\\$a" \
    --maf \${MAF_THRESH} \
    --hwe \${HWE_THRESH} \
    --exclude \${INFO_FAILED_LIST} \
    --extract "${VARIANT_LIST_MNT_PATH}" \
    --snps-only \
    --make-pgen \
    --out "\${OUTPUT_PREFIX_NEW}"
    
if [ -f "\${OUTPUT_PREFIX_NEW}.pvar" ]; then
    count=\$((\$(wc -l < "\${OUTPUT_PREFIX_NEW}.pvar") - 1))
    
    # Clean up
    rm *snplist *log "\${CORRECTED_FAM_PATH}" "\${TEMP_SMALL_BGEN}"
    
    echo "Success: Chr $CHR. Wrote \$count passing variants to \${OUTPUT_PREFIX_NEW}.pgen"
else
    echo "Error: Chr $CHR filtering failed."
    exit 1
fi

EOF


# --- 6. Execute the dx run command ---
dx run swiss-army-knife \
  -icmd="$SCRIPT_CMD" \
  --name="$JOB_NAME" \
  --instance-type "$instance_type" \
  --priority="$priority" \
  --destination="$DESTINATION" \
  --brief \
  -y