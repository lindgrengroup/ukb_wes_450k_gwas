#!/bin/bash

# --- 1. Check Input ---
if [[ -z "$1" ]]; then
    echo "Error: Missing gene symbol." >&2
    echo "Usage: ./get_ensgid.sh [GENE_SYMBOL]" >&2
    exit 1
fi

GENE_SYMBOL=$1
ENSEMBL_GRCh37_GFF="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ensembl/Homo_sapiens.GRCh37.87.gff3.gz"

# --- 2. Find ENSGID ---
# Redirect stderr to null to suppress potential gunzip/pipe errors on broken pipes (exit early)
ensgid=$( cat "${ENSEMBL_GRCh37_GFF}" 2>/dev/null | gunzip -c | awk -F'\t' -v TARGET="$GENE_SYMBOL" '
!/^#/ && $3 == "gene" {
    split($9, attrs, ";");
    gene_name = "";
    gene_id = "";
    
    for (i in attrs) {
        # Extract Name to match input
        if (attrs[i] ~ /^Name=/) {
            sub("Name=", "", attrs[i]);
            gene_name = attrs[i];
        }
        # Extract ID (format usually ID=gene:ENSG...)
        if (attrs[i] ~ /^ID=/) {
            sub("ID=", "", attrs[i]);
            sub("gene:", "", attrs[i]); # Remove "gene:" prefix if present
            gene_id = attrs[i];
        }
    }
    
    if (gene_name == TARGET) {
        print gene_id;
        exit;
    }
}')

# --- 3. Check and Print Result ---
if [[ -z "$ensgid" ]]; then
    echo "Error: ENSGID for gene $GENE_SYMBOL not found." >&2
    exit 1
fi

echo "$ensgid"
