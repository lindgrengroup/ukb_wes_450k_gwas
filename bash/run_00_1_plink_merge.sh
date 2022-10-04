#!/bin/bash

set -u # throws error if variables are undefined

export project=`dx pwd`
WD="/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas"

readonly script="00_1_plink_merge.sh"

script_local="${WD}/bash/${script}"
script_dnax="/saige_pipeline/scripts/${script}"

source "/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file "${script_local}" "${script_dnax}"

dx run swiss-army-knife \
	-iin="/saige_pipeline/scripts/${script}" \
	-icmd="bash ${script}" \
	--name="plink_merge_pruned" \
	--instance-type "mem1_ssd1_v2_x8" \
	--destination="/saige_pipeline/data/00_set_up" \
	-y