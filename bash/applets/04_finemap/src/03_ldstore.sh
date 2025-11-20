#!/bin/bash
# 02_saige_variant_test 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    mkdir plink_files
    mkdir -p out/{output,log}_file
    
    ## Download input files
    for suffix in {bed,bim,fam}; do
      dx download "${plink_bfile}.${suffix}" -o "plink_files/bfile.${suffix}"
    done
    # Remove "chr" prefix from chromosome column in bim file
    sed -i 's/^chr//g' "plink_files/bfile.bim"

    dx download "${sparse_grm}" -o "sparse_grm"
    dx download "${sparse_grm_samples}" -o "sparse_grm_samples"

    dx download "${model_file}" -o "model_file"
    dx download "${variance_ratios}" -o "variance_ratios"

    # Download saige-1.1.6.1.tar.gz docker image from wes_450k/ukbb_meta/docker/
    dx download file-GGzB25jJg8Jzf3gQGxY0Jy4X 
    docker load --input saige-1.1.6.1.tar.gz

    ## Run script
    docker_dir="/mnt/home"
    docker run \
      -e chrom="${chrom}"  \
      -e output_prefix="${output_prefix}"  \
      -v "${WD}":"${docker_dir}" \
      wzhou88/saige:1.1.6.1 /bin/bash -c """step2_SPAtests.R  \
        --bedFile="${docker_dir}/plink_files/bfile.bed" \
        --bimFile="${docker_dir}/plink_files/bfile.bim" \
        --famFile="${docker_dir}/plink_files/bfile.fam" \
        --AlleleOrder=ref-first \
        --chrom="${chrom}" \
        --minMAF=0  \
        --minMAC=20  \
        --GMMATmodelFile="${docker_dir}/model_file"  \
        --varianceRatioFile="${docker_dir}/variance_ratios"  \
        --LOCO=FALSE  \
        --is_Firth_beta=TRUE  \
        --pCutoffforFirth=0.1 \
        --is_output_moreDetails=TRUE  \
        --sparseGRMFile="${docker_dir}/sparse_grm"  \
        --sparseGRMSampleIDFile="${docker_dir}/sparse_grm_samples" \
        --is_fastTest=TRUE  \
        --SAIGEOutputFile="${docker_dir}/${output_prefix}.tsv" 2>&1 | tee "${docker_dir}/${output_prefix}.log"
        """

    # docker run \
    #   -e chrom="${chrom}"  \
    #   -e output_prefix="${output_prefix}"  \
    #   -v "${WD}":"${docker_dir}" \
    #   wzhou88/saige:1.1.6.1 /bin/bash -c """echo 'test run' 2>&1 | tee "${docker_dir}/${output_prefix}.log"; touch ${docker_dir}/${output_prefix}.tsv"""

    gzip "${output_prefix}.tsv" 
    mv "${output_prefix}.tsv.gz" "out/output_file/" 
    mv "${output_prefix}.log" "out/log_file/"

    dx-upload-all-outputs
}