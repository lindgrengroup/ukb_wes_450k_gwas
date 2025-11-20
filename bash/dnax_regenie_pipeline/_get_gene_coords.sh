#!/bin/bash

# --- 1. Check Input ---
if [[ -z "$1" ]]; then
    echo "Error: Missing gene symbol." >&2
    echo "Usage: ./get_gene_coords.sh [GENE_SYMBOL]" >&2
    exit 1
fi

GENE_SYMBOL=$1
#ENSEMBL_GRCh37_GFF="/resources/ensembl/Homo_sapiens.GRCh37.87.gff3.gz" # DNAnexus file path
ENSEMBL_GRCh37_GFF="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ensembl/Homo_sapiens.GRCh37.87.gff3.gz" # Local file

# --- 2. Find Gene Coordinates ---
# Redirect stderr to null so we don't see the Python panic
gene_info=$( cat "${ENSEMBL_GRCh37_GFF}" 2>/dev/null | gunzip -c | awk -F'\t' -v TARGET="$GENE_SYMBOL" '
!/^#/ && $3 == "gene" {
    split($9, attrs, ";");
    gene_name = "";
    for (i in attrs) {
        if (attrs[i] ~ /^Name=/) {
            sub("Name=", "", attrs[i]);
            gene_name = attrs[i];
        }
    }
    if (gene_name == TARGET) {
        print $1, $4, $5; # Print CHR START STOP
        exit;
    }
}')

# --- 3. Check and Print Result ---
if [[ -z "$gene_info" ]]; then
    echo "Error: Gene $GENE_SYMBOL not found or GFF file read error." >&2
    exit 1
fi

# Print the found coordinates to standard output
echo "$gene_info"
