#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="00_2_basic_qc.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}

## OPTIONS
#readonly anc="EUR" # Genetic ancestry group (options: "EUR", "AFR", "EAS", "SAS", "AMR")
for anc in {"EAS","SAS","AMR"}; do
  #readonly BFILE="ukb22418_b0_v2.autosomes" # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")
  
  echo "Running REGENIE basic QC on $anc"
  
  dx run swiss-army-knife \
  	-iin="${script_dnax}" \
  	-icmd="bash ${script} $anc" \
  	--name="regenie_basic_qc_${anc}" \
  	--instance-type "mem3_ssd1_v2_x8" \
  	--priority="high" \
  	--destination="/nbaya/regenie/data/genotypes/" \
  	--brief \
  	-y
  
done
