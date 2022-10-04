#!/bin/bash

export project=`dx pwd`

set -u # throws error if variables are undefined

readonly script="00_1_plink_merge.sh"

dx run swiss-army-knife \
	-iin="/saige_pipeline/scripts/${script}" \
	-icmd="bash ${script}" \
	--name="plink merge pruned" \
	--instance-type "mem1_ssd1_v2_x8" \
	--destination="/saige_pipeline/data/00_sparse_grm" \
	-y