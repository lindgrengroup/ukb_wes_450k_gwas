#!/bin/bash

set -u # throws error if variables are undefined

export project=`dx pwd`

readonly script="00_2_create_sparse_grm.sh"

WD="/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas"
script_local="${WD}/bash/${script}"
script_dnax="/saige_pipeline/scripts/${script}"

source "/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file "${script_local}" "${script_dnax}"

bfile="/saige_pipeline/data/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned"

dx run swiss-army-knife \
	--name="00_2_create_sparse_grm" \
	-iin="/saige_pipeline/scripts/${script}" \
	-iin="${bfile}.bed" \
	-iin="${bfile}.bim" \
	-iin="${bfile}.fam" \
	-icmd="bash ${script}" \
	-iimage_file="/saige_pipeline/docker_images/saige_1.0.9.tar.gz" \
	--instance-type "mem1_ssd1_v2_x72" \
	--destination="/saige_pipeline/data/00_set_up" \
	--brief \
	-y 
