#!/bin/bash

set -u # throws error if variables are undefined

#WD="/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas"
WD="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"
source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh" # needed for `upload_file` function

# Destination of output files
outdir="/saige_pipeline/data/sentinel_variants/wes_merged_w_sentinels-max_log10bf-include_ties/"

readonly script="00_merge_wes_with_sentinel_variants.sh"

script_local="${WD}/bash/${script}"
script_dnax="/saige_pipeline/scripts/${script}"

# Upload script to DNAnexus
upload_file "${script_local}" "${script_dnax}"

for chrom in {1..20}; do
  echo "Running merge_wes_with_sentinel for chr${chrom}"
  
  dx run swiss-army-knife \
  	-iin="/saige_pipeline/scripts/${script}" \
  	-icmd="bash ${script} ${chrom}" \
  	--name="merge_wes_with_sentinel-c${chrom}" \
  	--instance-type "mem2_ssd2_v2_x16" \
  	--priority="high" \
  	--destination="${outdir}" \
  	--brief \
  	-y

done
