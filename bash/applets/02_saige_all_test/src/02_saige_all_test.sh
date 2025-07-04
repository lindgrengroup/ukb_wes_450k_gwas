#!/bin/bash
# 02_saige_variant_test 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    mkdir -p out/{gene_assoc,marker_list,variant_assoc,log_file}
    
    dx-mount-all-inputs --except bed --except bim
    dx download "$bed" -o "${WD}/bfile.bed"
    dx download "$bim" -o "${WD}/bfile.bim"
    
    if [ $( head -1 bfile.bim | cut -f1 ) = "X" ]; then 
      # If chrX, need to fill in missing varids
      awk '{
        
        varid="chr"$1":"$4":"$6":"$5  

        print $1,varid,$3,$4,$5,$6

      }' "${WD}/bfile.bim" > "${WD}/bfile.bim-tmp"

      mv "${WD}/bfile.bim-tmp" "${WD}/bfile.bim"
    else
      # Remove "chr" prefix from chromosome column in bim file
      sed -i 's/^chr//g' "${WD}/bfile.bim"
    fi

    dx download "${group_file}" -o "group_file.gz"
    gunzip -c group_file.gz > "${WD}/group_file.txt"

    # NOTE: Started using v1.1.9 on 2023-07-18. SAIGE v1.1.8 is the first version that works correctly for case-control traits.
    # Previously I was using v1.1.6.3 (and before that, v1.1.6.1)
    docker pull wzhou88/saige:1.1.9

    ## Run script
    docker run \
      -e HOME=${WD}  \
      -v ${WD}/:$HOME/ \
      wzhou88/saige:1.1.9 \
      step2_SPAtests.R  \
        --bedFile ${HOME}/bfile.bed \
        --bimFile ${HOME}/bfile.bim \
        --famFile ${HOME}/in/fam/* \
        --LOCO FALSE  \
        --AlleleOrder ref-first \
        --minMAF 0  \
        --minMAC 0.5  \
        --GMMATmodelFile ${HOME}/in/model_file/*  \
        --varianceRatioFile ${HOME}/in/variance_ratios/*  \
        --sparseGRMFile ${HOME}/in/sparse_grm/*  \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm_samples/* \
        --groupFile ${HOME}/group_file.txt \
        --annotation_in_groupTest "pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous" \
        --maxMAF_in_groupTest 0.0001,0.001,0.01 \
        --is_output_markerList_in_groupTest TRUE \
        --is_single_in_groupTest TRUE \
        --is_Firth_beta TRUE    \
        --pCutoffforFirth 0.05 \
        --is_fastTest TRUE  \
        --SAIGEOutputFile ${HOME}/${output_prefix}.tsv \
        2>&1 | tee ${HOME}/${output_prefix}.log

    gzip "${output_prefix}.tsv" 
    mv "${output_prefix}.tsv.gz" "out/gene_assoc/" 
    gzip "${output_prefix}.tsv.markerList.txt" 
    mv "${output_prefix}.tsv.markerList.txt.gz" "out/marker_list/"
    gzip "${output_prefix}.tsv.singleAssoc.txt" 
    mv "${output_prefix}.tsv.singleAssoc.txt.gz" "out/variant_assoc/"
    mv "${output_prefix}.log" "out/log_file/"

    dx-upload-all-outputs
}
