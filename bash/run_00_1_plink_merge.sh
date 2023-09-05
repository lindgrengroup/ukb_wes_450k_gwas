#!/bin/bash

set -u # throws error if variables are undefined

WD="/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas"

readonly script="00_1_plink_merge.sh"

script_local="${WD}/bash/${script}"
script_dnax="/saige_pipeline/scripts/${script}"

source "/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file "${script_local}" "${script_dnax}"

## OPTIONS
# dataset_to_merge:
# - 'pruned'
# - 'for_vr'
readonly dataset_to_merge="for_vr"

# pop:
# - 'allpop'
# - 'eur'
readonly pop='allpop'

echo "Running PLINK merge on: ${pop}-${dataset_to_merge}"

dx run swiss-army-knife \
	-iin="/saige_pipeline/scripts/${script}" \
	-icmd="bash ${script} ${pop} ${dataset_to_merge}" \
	--name="plink_merge_${pop}_${dataset_to_merge}" \
	--instance-type "mem1_ssd1_v2_x8" \
	--priority="low" \
	--destination="/saige_pipeline/data/00_set_up" \
	--brief \
	-y