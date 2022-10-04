#!/bin/bash

export project=`dx pwd`

set -u # throws error if variables are undefined

readonly script="00_2_create_sparse_grm.sh"

# instance_type="mem1_ssd1_v2_x36"
instance_type="mem1_ssd1_v2_x72"

dx run swiss-army-knife \
	--name="create sparse grm" \
	-iin="/saige_pipeline/scripts/${script}" \
	-iin="/saige_pipeline/data/00_sparse_grm/ukb_array_pruned.bed" \
	-iin="/saige_pipeline/data/00_sparse_grm/ukb_array_pruned.bim" \
	-iin="/saige_pipeline/data/00_sparse_grm/ukb_array_pruned.fam" \
	-icmd="bash ${script}" \
	-iimage_file="/docker_images/saige_1.0.9.tar.gz" \
	--instance-type "${instance_type}" \
	--destination="/saige_pipeline/data/00_sparse_grm" \
	-y 
