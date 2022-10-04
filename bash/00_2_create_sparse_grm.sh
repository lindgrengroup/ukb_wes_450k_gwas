#!/bin/bash 
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

bfile="ukb_array_pruned"
out="test-ukb_array"

createSparseGRM.R       \
    --plinkFile="${bfile}" \
    --nThreads=72  \
    --outputPrefix="${out}"       \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05