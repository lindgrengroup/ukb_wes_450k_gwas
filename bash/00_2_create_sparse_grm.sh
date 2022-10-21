#!/bin/bash 
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

set -e # Stop job if any command fails

bfile="ukb_array.wes_450k_qc_pass_eur.pruned"
out="ukb_array.wes_450k_qc_pass_eur.pruned"

createSparseGRM.R       \
    --plinkFile="${bfile}" \
    --nThreads=72  \
    --outputPrefix="${out}"       \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05