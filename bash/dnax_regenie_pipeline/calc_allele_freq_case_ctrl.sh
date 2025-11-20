#!/bin/bash

# --- Script Configuration ---
# Exit immediately if a command exits with a non-zero status.
set -e

# --- 1. Validate Input ---
# Check if the phenotype group argument is provided.
if [ "$#" -ne 4 ]; then
  echo "Error: You must provide a phenotype group name, phenotype, chromosome and output prefix."
  echo "Usage: $0 <phenotype_group> <phenotype> <chrom> <output>"
  exit 1
fi

PHENOTYPE_GROUP=$1
PHENO=$2
CHROM=$3
OUT=$4

# --- 2. Define Paths and Variables ---
# Centralize file paths for easier management.
BASE_DIR="/mnt/project"
BFILE_DIR="${BASE_DIR}/Barney/wes/sample_filtered"
PHENO_DIR="${BASE_DIR}/saige_pipeline/data/phenotypes"

# Input files
PHENO_TSV_PATH="${PHENO_DIR}/ukb_wes.phenos_and_covariates.${PHENOTYPE_GROUP}.tsv.gz"
PHENO_LIST_PATH="${PHENO_DIR}/phenotype_list.${PHENOTYPE_GROUP}.txt"
BASE_FAM_TEMPLATE="${BFILE_DIR}/ukb_wes_450k.qced.chr21.fam" # Template from any chromosome is fine

# Output directory
OUTPUT_DIR="./"

# --- 3. Setup Environment ---
# Create the output directory. The -p flag prevents errors if it already exists.
mkdir -p "$OUTPUT_DIR"

# Download and unzip PLINK if the executable is not found in the current directory.
if [ ! -f ./plink ]; then
  echo "PLINK not found. Downloading..."
  wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip
  unzip -q *zip plink
  rm *zip
  echo "PLINK downloaded successfully."
fi

# Check if the required phenotype files exist before starting.
if [ ! -f "$PHENO_TSV_PATH" ] || [ ! -f "$PHENO_LIST_PATH" ]; then
  echo "Error: Phenotype data not found at:"
  echo "TSV: ${PHENO_TSV_PATH}"
  echo "List: ${PHENO_LIST_PATH}"
  exit 1
fi

# --- 4. Create Helper Script to Generate FAM Files ---
# A Python script is better for this data manipulation task.
# We create it on-the-fly using a 'here document'.
cat << EOF > generate_fam.py
import pandas as pd
import sys

# Arguments from bash script
pheno_name = sys.argv[1]
pheno_tsv_path = sys.argv[2]
base_fam_path = sys.argv[3]
output_fam_path = f"tmp-{pheno_name}.fam"

try:
    # 1. Load data
    pheno_df = pd.read_csv(pheno_tsv_path, sep='\\t', usecols=['IID', pheno_name])
    base_fam = pd.read_csv(base_fam_path, sep=r'\\s+', header=None,
                         names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])

    # 2. Merge phenotype data with the base FAM file
    # Use a left merge to keep all individuals from the original FAM file
    merged_fam = base_fam.drop(columns=['PHENO']).merge(
        pheno_df, on='IID', how='left'
    )

    # 3. Recode phenotype for PLINK
    # Original: 0=control, 1=case. PLINK standard: 1=control, 2=case.
    merged_fam[pheno_name] += 1

    print(f'{pheno_name} (n={merged_fam[pheno_name].count()})')
    print(f'* cases: {(merged_fam[pheno_name]==2).sum()}') # https://www.cog-genomics.org/plink/1.9/formats#fam
    print(f'* ctrls: {(merged_fam[pheno_name]==1).sum()}')

    # Missing phenotypes should be -9 for PLINK to ignore them
    merged_fam[pheno_name] = merged_fam[pheno_name].fillna(-9)

    # Ensure the phenotype column is integer type
    merged_fam[pheno_name] = merged_fam[pheno_name].astype(int)

    # 4. Save the new FAM file
    # Ensure columns are in the correct order for the .fam format
    final_fam = merged_fam[['FID', 'IID', 'PAT', 'MAT', 'SEX', pheno_name]]
    final_fam.to_csv(output_fam_path, sep=' ', header=False, index=False)

    # print(f"Successfully generated: {output_fam_path}")

except Exception as e:
    print(f"Error processing phenotype '{pheno_name}': {e}", file=sys.stderr)
    sys.exit(1)
EOF

# --- 5. Main Processing Loop ---
echo "✅ Starting analysis for phenotype group=${PHENOTYPE_GROUP} ${PHENO} chr${CHROM}"

# Read the phenotype list file line by line

# Generate the temporary .fam file needed for this phenotype
python3 generate_fam.py "${PHENO}" "$PHENO_TSV_PATH" "$BASE_FAM_TEMPLATE"
TEMP_FAM_FILE="tmp-${PHENO}.fam"

# Use local 
BFILE_PREFIX="/mnt/project/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr${CHROM}"

echo "Using PLINK files: ${BFILE_PREFIX}.{bed,bim}"

# Check that the required chromosome files exist before running PLINK
if [ ! -f "${BFILE_PREFIX}.bed" ] || [ ! -f "${BFILE_PREFIX}.bim" ]; then
  echo "    ⚠️ Warning: BFILE for chr${chrom} (${BFILE_PREFIX}.{bed,bim}) not found. Skipping."
  continue
fi

# Execute PLINK command
# Note the `\` at the end of each line allows the command to be split for readability.
# Include '--allow-no-sex' to include individuals with missing sex
./plink \
  --bed "${BFILE_PREFIX}.bed" \
  --bim "${BFILE_PREFIX}.bim" \
  --fam "${TEMP_FAM_FILE}" \
  --freq case-control \
  --allow-no-sex \
  --out "${OUT}" > /dev/null # Suppress verbose PLINK output

# Clean up the temporary FAM file for the current phenotype
rm "${TEMP_FAM_FILE}"


# --- 6. Final Cleanup ---
rm generate_fam.py plink *nosex *log
echo "🎉 Analysis complete: ${out}."
