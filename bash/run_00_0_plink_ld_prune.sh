#!/bin/bash

export project=`dx pwd`

set -u # throws error if variables are undefined

readonly chrom=$1
readonly script="00_0_plink_ld_prune.sh"

dx run swiss-army-knife \
	-iin="/saige_pipeline/scripts/${script}" \
	-icmd="bash ${script} ${chrom}" \
	--name="plink ld prune c${chrom}" \
	--instance-type "mem1_ssd1_v2_x2" \
	--destination="/saige_pipeline/data/00_sparse_grm" \
	-y