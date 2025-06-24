#!/bin/bash
# 01_saige_fit_null 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    # mkdir -p "plink_files"
    mkdir -p out/{model_file,variance_ratios}

    dx-mount-all-inputs
    
    # Download phenotype file
    zcat ${HOME}/in/pheno_file/* > "pheno_file"

    # NOTE: Started using v1.1.9 on 2023-07-18. SAIGE v1.1.8 is the first version that works correctly for case-control traits. 
    # Previously I was using v1.1.6.3 (and before that, v1.1.6.1)
    docker pull wzhou88/saige:1.1.9 

    # Get number of threads
    n_threads=$(( $(nproc --all) - 1 ))

    if [ ${sex} == "female" ] || [ ${sex} == "male" ]; then
      sex_specific_sample_list="sample_ids.txt"
      dx download "project-GBvkP10Jg8Jpy18FPjPByv29:/saige_pipeline/data/covariates/ukb_sex.tsv.gz"
      if [ ${sex} == "female" ]; then
        is_female=1
      else
        is_female=0
      fi

      zcat ukb_sex.tsv.gz \
      | awk -v is_female="${is_female}" '$2==is_female { print $1 }' > ${sex_specific_sample_list}
      sample_id_include_flag="--SampleIDIncludeFile=${WD}/${sex_specific_sample_list}"
    else
      sample_id_include_flag=""
    fi

    head ${sex_specific_sample_list}

    # Get inverse-normalize flag if trait_type=="quantitative"
    if [ ${trait_type} == "quantitative" ]; then
      if [ ${use_irnt} == "true" ]; then
        use_irnt_flag="--invNormalize=TRUE"
      else
        use_irnt_flag="--invNormalize=FALSE"
      fi
    fi

    ## Run script
    docker run \
      -e HOME=${WD}  \
      -e pheno_col="${pheno_col}" \
      -e trait_type="${trait_type}" \
      -e inv_normalize_flag=${inv_normalize_flag} \
      -e output_prefix="${output_prefix}"  \
      -e n_threads="${n_threads}" \
      -v ${WD}/:$HOME/ \
      wzhou88/saige:1.1.9 step1_fitNULLGLMM.R  \
        --bedFile ${HOME}/in/plink_for_vr_bed/* \
        --bimFile ${HOME}/in/plink_for_vr_bim/* \
        --famFile ${HOME}/in/plink_for_vr_fam/* \
        --sparseGRMFile ${HOME}/in/sparse_grm/* \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm_samples/*  \
        --useSparseGRMtoFitNULL=TRUE \
        --phenoFile ${HOME}/pheno_file \
        ${sample_id_include_flag} \
        --skipVarianceRatioEstimation FALSE \
        --phenoCol "${pheno_col}" \
        --covarColList "${covar_col_list}" \
        --qCovarColList="${qcovar_col_list}"  \
        --sampleIDColinphenoFile="IID" \
        --isCateVarianceRatio=FALSE \
        --traitType=${trait_type} \
        ${use_irnt_flag} \
        --outputPrefix="${HOME}/${output_prefix}" \
        --IsOverwriteVarianceRatioFile=TRUE \
        --nThreads=${n_threads} \
        --LOCO=TRUE 
        #--isLowMemLOCO=TRUE
        # Added isLowMemLOCO=TRUE on 10 Oct 2024 (NOTE: --LOCO was already TRUE before this addition)

    mv *.rda out/model_file/
    mv *.varianceRatio.txt out/variance_ratios/

    dx-upload-all-outputs
}

