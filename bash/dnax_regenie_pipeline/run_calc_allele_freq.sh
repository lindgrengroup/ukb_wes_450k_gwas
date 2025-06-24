#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="calc_allele_freq.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}

## OPTIONS
anc="SAS" # Genetic ancestry group (options: "EUR", "AFR", "EAS", "SAS", "AMR")
#for anc in {"EAS","SAS","AMR"}; do

PHENOCOL="Height"

DESTINATION_DIR="/nbaya/regenie/data/allele_freq" # DNAnexus destination directory

for chrom in {7..23}; do
  if [ $chrom -eq 23 ]; then
    chrom="X"
  fi
  #readonly BFILE="ukb22418_b0_v2.autosomes" # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")

  OUTPUT_FILE="ukb_wes_450k.qced.chr${chrom}.${anc}.${PHENOCOL}.frq"
  FULL_DX_PATH="${DESTINATION_DIR}/${OUTPUT_FILE}"

  if ! dx ls "${FULL_DX_PATH}" &> /dev/null; then
  
    echo "Running REGENIE basic QC on $anc for chr$chrom"
    
    dx run swiss-army-knife \
    	-iin="${script_dnax}" \
    	-icmd="bash ${script} ${anc} ${PHENOCOL} ${chrom}" \
    	--name="calc_allele_freq_${anc}_${PHENOCOL}_chr${chrom}" \
    	--instance-type "mem2_ssd1_v2_x8" \
    	--priority="high" \
    	--destination="/nbaya/regenie/data/allele_freq" \
    	--brief \
    	-y
  else
    echo "File ${OUTPUT_FILE} already exists. Skipping chromosome ${chrom}."
  fi

  
done
