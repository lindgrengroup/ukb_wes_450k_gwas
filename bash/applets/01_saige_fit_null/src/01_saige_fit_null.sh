#!/bin/bash
# 01_saige_fit_null 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    mkdir -p "plink_files"
    mkdir -p out/{model_file,variance_ratios}
    
    ## Download input files
    for suffix in {bed,bim,fam}; do
      dx download "${plink_for_vr_bfile}.${suffix}" -o "plink_files/bfile.${suffix}"
    done
    dx download "$sparse_grm" -o "sparse_grm"
    dx download "$sparse_grm_samples" -o "sparse_grm_samples"
    
    dx download "$pheno_file" -o "pheno_file.tsv.gz"
    zcat "pheno_file.tsv.gz" > pheno_file

    # Download saige-1.1.6.1.tar.gz docker image from wes_450k/ukbb_meta/docker/
    dx download file-GGzB25jJg8Jzf3gQGxY0Jy4X 
    docker load --input saige-1.1.6.1.tar.gz

    ## Run script
    docker_dir="/mnt/home"
    docker run \
      -e pheno_col="${pheno_col}" \
      -e output_prefix="${output_prefix}"  \
      -v "${WD}":"${docker_dir}" \
      wzhou88/saige:1.1.6.1 step1_fitNULLGLMM.R  \
        --sparseGRMFile="${docker_dir}/sparse_grm" \
        --sparseGRMSampleIDFile="${docker_dir}/sparse_grm_samples"  \
        --useSparseGRMtoFitNULL=TRUE  \
        --plinkFile="${docker_dir}/plink_files/bfile" \
        --phenoFile="${docker_dir}/pheno_file" \
        --skipVarianceRatioEstimation=FALSE \
        --phenoCol="${pheno_col}" \
        --covarColList="age,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15,pc16,pc17,pc18,pc19,pc20,pc21,age2,is_female_age,is_female_age2,is_female,assessment_centre,genotyping_batch" \
        --qCovarColList="is_female,assessment_centre,genotyping_batch"  \
        --sampleIDColinphenoFile="IID" \
        --traitType="quantitative"        \
        --outputPrefix="${docker_dir}/${output_prefix}" \
        --IsOverwriteVarianceRatioFile=TRUE \
        --invNormalize=TRUE

    mv *.rda out/model_file/
    mv *.varianceRatio.txt out/variance_ratios

    dx-upload-all-outputs
}