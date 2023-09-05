#!/bin/bash
# 02_saige_set_test 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    mkdir -p out/{output,log}_file
    
    dx-mount-all-inputs --except bed --except bim
    dx download "$bed" -o "${WD}/bfile.bed"
    dx download "$bim" -o "${WD}/bfile.bim"

    if [ $( head -1 bfile.bim | cut -f1 ) = "chrX" ]; then 
      # If chrX, need to fill in missing varids
      awk '{
        
        varid=$1":"$4":"$6":"$5  

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
      wzhou88/saige:1.1.6.3 \
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
        --groupFile ${HOME}/in/group_file/* \
        --annotation_in_groupTest "pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous" \
        --maxMAF_in_groupTest 0.0001,0.001,0.01 \
        --is_output_markerList_in_groupTest TRUE \
        --is_fastTest TRUE \
        --SAIGEOutputFile ${HOME}/${output_prefix}.tsv 2>&1 | tee ${HOME}/${output_prefix}.log
        
    # --groupFile ${HOME}/group_file.txt \
    # pLoF,damaging_missense,synonymous,non_coding,pLoF:damaging_missense,pLoF:damaging_missense:synonymous,pLoF:damaging_missense:synonymous:non_coding
    # pLoF,damaging_missense,other_missense,synonymous,non_coding,pLoF:damaging_missense,pLoF:damaging_missense:other_missense,pLoF:damaging_missense:other_missense:synonymous,pLoF:damaging_missense:other_missense:synonymous:non_coding
    # pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense,pLoF:damaging_missense:other_missense:synonymous
    # pLoF,damaging_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:synonymous

    gzip "${output_prefix}.tsv" 
    mv "${output_prefix}.tsv.gz" "out/output_file/" 
    mv "${output_prefix}.log" "out/log_file/"

    dx-upload-all-outputs
}

