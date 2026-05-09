#!/bin/bash

# --- 1. Local Setup & Lookups ---

# Your input grouped gene ID
# Options:
# - "t2d-high_outlier_canonical_genes"
# - "cad-low_outlier_canonical_genes"
# - "osteoporosis-high_outlier_canonical_genes"
GROUP_ID="$1"

# Define gene array
if [[ ${GROUP_ID} == "t2d-high_outlier_canonical_genes" ]]; then
	GENE_ARRAY=( HNF1A HNF1B HNF4A GCK )
elif [[ ${GROUP_ID} == "cad-low_outlier_canonical_genes" ]]; then
	GENE_ARRAY=( APOB PCSK9 ANGPTL3 )
elif [[ ${GROUP_ID} == "osteoporosis-high_outlier_canonical_genes" ]]; then
	GENE_ARRAY=( COL1A1 COL1A2 )
else
	echo "Invalid ${GROUP_ID}!" && exit 1
fi

# Initialize arrays to hold the looked-up data
ENSGID_LIST=()
CHROM_LIST=()

echo "Performing local lookups..."

for gene in "${GENE_ARRAY[@]}"; do
	# 1. Get ENSGID using your local function
	gid=$(./_get_gene_ensgid.sh "$gene" "38") # Use build 38 because UKB WES uses GRCh38

	# 2. Get Coords (Chrom) using your local function
	# Assuming _get_gene_coords returns: "9 12345 67890" (chrom start end)
	# We capture the output into an array to grab the first element
	coords_output=($(./_get_gene_coords.sh "$gene"))
	chrom="${coords_output[0]}"

	# Check if lookup failed (optional safety)
	if [[ -z "$gid" || -z "$chrom" ]]; then
		echo "Error: Lookup failed for $gene"
		exit 1
	fi

	# Append to lists
	ENSGID_LIST+=("$gid")
	CHROM_LIST+=("$chrom")
done

# Convert arrays to space-separated strings for injection into SAK
# These variables will be "burned" into the script below
STR_GENES="${GENE_ARRAY[*]}"
STR_IDS="${ENSGID_LIST[*]}"
STR_CHROMS="${CHROM_LIST[*]}"

echo "Lookups complete."
echo "Genes:  $STR_GENES"
echo "IDs:    $STR_IDS"
echo "Chroms: $STR_CHROMS"

# --- 2. Define the SAK Job Command ---

# We inject the strings created in Cell 1 directly into the variables below
read -r -d '' SCRIPT_CMD <<EOF
set -e

# --- Part A: Reconstruct Arrays on Worker ---
# We transform the space-separated strings back into Bash arrays
# Note: These variables are from the OUTER script, so they are NOT escaped
input_genes=($STR_GENES)
input_ids=($STR_IDS)
input_chroms=($STR_CHROMS)

# Define directories
DATA_DIR="/mnt/project/Barney/wes/sample_filtered"
ANNOT_DIR="/mnt/project/nbaya/regenie/data/annotations/v6"
OUT_PREFIX="${GROUP_ID}"

mkdir -p out/temp
MERGE_LIST="out/merge_list.txt"
touch \${MERGE_LIST}

echo "Starting processing for \${#input_genes[@]} genes..."

# --- Part B: Parallel Array Loop ---
# We loop using index 'i' to access gene, id, and chrom simultaneously
for ((i=0; i<\${#input_genes[@]}; i++)); do

		# Access array elements using the index
		gene="\${input_genes[\$i]}"
		ensgid="\${input_ids[\$i]}"
		chrom="\${input_chroms[\$i]}" # This is the specific chrom for this gene

		echo "Processing \$gene (\$ensgid) on Chr \${chrom}..."

		# 1. Define Paths
		# We now know exactly which annotation file to use based on the chrom passed in
		annot_file="\${ANNOT_DIR}/regenie_annotations.v6.chr\${chrom}.txt"
		var_list="out/temp/\${gene}_\${ensgid}.variants"

		# 2. Extract Variants
		# We grep the specific ENSGID from the specific chromosome file
		grep "\$ensgid" "\${annot_file}" | awk '{print \$1}' > \${var_list}

		# 3. PLINK Extract
		bfile_path="\${DATA_DIR}/ukb_wes_450k.qced.chr\${chrom}"
		gene_out="out/temp/\${gene}_chr\${chrom}"

		# Check for bfile existence
		if [[ ! -f "\${bfile_path}.bed" ]]; then
				echo "  Error: Bfile not found: \${bfile_path}"
				continue
			fi

		plink \
				--bfile "\${bfile_path}" \
				--extract "\${var_list}" \
				--make-bed \
				--out "\${gene_out}" \

		# 4. Add to Merge List
		echo "\${gene_out}" >> \${MERGE_LIST}

	done

# --- Part C: Merge ---
# Only run merge if we successfully added files to the list
if [[ -s \${MERGE_LIST} ]]; then
		echo "Merging datasets..."
		plink \
				--merge-list \${MERGE_LIST} \
				--make-bed \
				--out "out/\${OUT_PREFIX}_merged"

		# Move final results to upload area
		mv out/\${OUT_PREFIX}_merged.bed .
		mv out/\${OUT_PREFIX}_merged.bim .
		mv out/\${OUT_PREFIX}_merged.fam .
		mv out/\${OUT_PREFIX}_merged.log .
	else
		echo "Merge list empty. No files generated."
	fi

# Cleanup
rm -rf out/
echo "Job Complete."
EOF

# --- 3. Run the SAK Job ---
instance_type="mem1_ssd1_v2_x8"
destination="nbaya/outliers/data/ukb_wes_grouped_genes"
dx mkdir -p $destination

dx run swiss-army-knife \
	-icmd="$SCRIPT_CMD" \
	--name="get_grouped_gene_wes450k_bfiles_${GROUP_ID}" \
	--instance-type "$instance_type" \
	--destination="$destination" \
	--brief \
  --priority="high" \
	-y
