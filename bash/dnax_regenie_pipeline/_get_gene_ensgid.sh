#!/bin/bash

# --- 1. Check Input ---
if [[ -z "$1" ]]; then
    echo "Error: Missing gene symbol." >&2
    echo "Usage: ./get_ensgid.sh [GENE_SYMBOL]" >&2
    exit 1
fi

GENE_SYMBOL=$1
DEFAULT_GRCH_VERSION=$2 # Options: 37, 38

ENSEMBL_GRCh37_GFF="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ensembl/Homo_sapiens.GRCh37.87.gff3.gz"
ENSEMBL_GRCh38_GFF="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ensembl/Homo_sapiens.GRCh38.115.gff3.gz"

# --- 2. Find ENSGID ---
# Redirect stderr to null to suppress potential gunzip/pipe errors on broken pipes (exit early)
get_ensgid() {
  _ensembl_gff=$1
  ensgid=$( cat "${_ensembl_gff}" 2>/dev/null | gunzip -c | awk -F'\t' -v TARGET="$GENE_SYMBOL" '
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

  echo $ensgid
}


ensgid_grch37=$( get_ensgid ${ENSEMBL_GRCh37_GFF} )
ensgid_grch38=$( get_ensgid ${ENSEMBL_GRCh38_GFF} )

if [[ "${ensgid_grch37}" != "${ensgid_grch38}" ]]; then
  echo "Warning: Ensembl gene ID is not the same between GRCh37 ($ensgid_grch37) and GRCh38 ($ensgid_grch38) for ${GENE_SYMBOL}" >&2

  if [[ ! -n $DEFAULT_GRCH_VERSION ]]; then # If no default set
    echo "To set a default GRCh version ('37' or '38') to use, include it as a second arg (./_get_gene_ensgid.sh <GENE> <DEFAULT_GRCH_VERSION>). Exiting." >&2
    exit 1
  elif [[ "${DEFAULT_GRCH_VERSION}" == "37" ]]; then
    echo "Defaulting to GRCh${DEFAULT_GRCH_VERSION} version of ENSGID" >&2
    ensgid_final="$ensgid_grch37"
  elif [[ "${DEFAULT_GRCH_VERSION}" == "38" ]]; then
    echo "Defaulting to GRCh${DEFAULT_GRCH_VERSION} version of ENSGID" >&2
    ensgid_final="$ensgid_grch38"
  fi
else
  ensgid_final=$ensgid_grch37 # If same, use GRCh37 (arbitrarily chosen)
fi

echo "$ensgid_final"
