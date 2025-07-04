#!/bin/bash
# 02_saige_all_test_multipheno 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    # mkdir -p out/{output,log}_files
    mkdir -p out/{gene_assoc,marker_list,variant_assoc,log_file}

    # # dx_wd=$( project-GBvkP10Jg8Jpy18FPjPByv29:/ )
    saige_data_dir="project-GBvkP10Jg8Jpy18FPjPByv29://saige_pipeline/data"
    output_dir="${saige_data_dir}/02_saige_set_test${catevr}/${phenotype_group}/${pop}/${sex}"
    
    # [INPUT] SAIGE null model file directory
    if [ ${use_irnt} == "true" ]; then
      model_suffix=""
    else
      model_suffix="-no_irnt"
    fi
    fit_null_output_prefix="${saige_data_dir}/01_fit_null${model_suffix}${catevr}/${phenotype_group}/${pop}/${sex}"

    # [INPUT] Pheno list
    if [[ "${phenotype_group}" == *"qced_biomarkers"* ]]; then
      # Always use the same list of biomarkers, with no suffixes related to tail status
      phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.qced_biomarkers.txt" ) )
    else
      phenos=( $( dx cat "${saige_data_dir}/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
    fi
    # # n_phenos=${#phenos[@]} # total number of phenotypes in list


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

    pheno_idx_start=$((pheno_idx_start))
    pheno_idx_stop=$((pheno_idx_stop))
    
    for pheno_idx in `seq $pheno_idx_start $pheno_idx_stop`; do

      pheno_col="${phenos[$pheno_idx]}" 
      
			# Only include IRNT-ed residualized phenotypes (specific to Mahalanobis phenotypes)
    	if [[ $phenotype_group == "mahalanobis_v2"* ]] && [[ $pheno_col != *"irnt"* ]]; then
    	  continue
    	fi

      gwas_id="${pheno_col}-${pop}-${sex}"

      # TEMPORARY (for testing LOCO option)
      #gwas_id="${gwas_id}-loco"
      loco_flag="--LOCO FALSE" # Default
      #loco_flag="--LOCO TRUE --chrom $chrom" # For testing LOCO option


      # [OUTPUT] Output files
      output_prefix="saige_all_test${catevr}.${gwas_id}.chr${chrom}"
      #output_prefix="not_fasttest-saige_all_test${catevr}.${gwas_id}.chr${chrom}"
      results_file="${output_dir}/${output_prefix}.tsv.gz"

      if [ $( dx ls -l ${results_file} 2> /dev/null | wc -l  ) -eq 0 ]; then
  
        echo "Starting $gwas_id"
        
        # # [INPUT] SAIGE null model files
        model_file="${fit_null_output_prefix}/${gwas_id}.rda"
        
        # TEMPORARY (for testing LOCO option)
        #model_file="${fit_null_output_prefix}/${gwas_id}_noLOCO.rda"


        variance_ratios="${fit_null_output_prefix}/${gwas_id}.varianceRatio.txt"
        dx download --overwrite "${model_file}" -o "${WD}/model_file.rda"
        dx download --overwrite "${variance_ratios}" -o "${WD}/variance_ratios.txt"

        # ## Run script
        docker run \
          -e HOME=${WD}  \
          -v ${WD}/:$HOME/ \
          wzhou88/saige:1.1.9 \
          step2_SPAtests.R  \
            --bedFile ${HOME}/bfile.bed \
            --bimFile ${HOME}/bfile.bim \
            --famFile ${HOME}/in/fam/* \
            ${loco_flag}  \
            --AlleleOrder ref-first \
            --minMAF 0  \
            --minMAC 0.5  \
            --GMMATmodelFile ${HOME}/model_file.rda  \
            --varianceRatioFile ${HOME}/variance_ratios.txt  \
            --sparseGRMFile ${HOME}/in/sparse_grm/*  \
            --sparseGRMSampleIDFile ${HOME}/in/sparse_grm_samples/* \
            --groupFile ${HOME}/in/group_file/* \
            --annotation_in_groupTest "pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous" \
            --maxMAF_in_groupTest 0.0001,0.001,0.01 \
            --is_output_markerList_in_groupTest TRUE \
            --is_fastTest TRUE \
            --SAIGEOutputFile ${HOME}/${output_prefix}.tsv 2>&1 | tee ${HOME}/${output_prefix}.log
      
        gzip "${output_prefix}.tsv" 
        mv "${output_prefix}.tsv.gz" "out/gene_assoc/" 
        gzip "${output_prefix}.tsv.markerList.txt" 
        mv "${output_prefix}.tsv.markerList.txt.gz" "out/marker_list/"
        gzip "${output_prefix}.tsv.singleAssoc.txt" 
        mv "${output_prefix}.tsv.singleAssoc.txt.gz" "out/variant_assoc/"
        mv "${output_prefix}.log" "out/log_file/"

        # gzip "${output_prefix}.tsv" 
        # mv "${output_prefix}.tsv.gz" "out/output_files/" 
        # mv "${output_prefix}.log" "out/log_files/"

        # dx upload "${output_prefix}.tsv.gz" --path ${results_file}
        # dx upload "${output_prefix}.log" --path "${output_dir}/"

        # Flag used for getting burden beta on the same scale as GWAS beta: --is_no_weight_in_groupTest TRUE \
      fi

    done

    # --groupFile ${HOME}/group_file.txt \
    # pLoF,damaging_missense,synonymous,non_coding,pLoF:damaging_missense,pLoF:damaging_missense:synonymous,pLoF:damaging_missense:synonymous:non_coding
    # pLoF,damaging_missense,other_missense,synonymous,non_coding,pLoF:damaging_missense,pLoF:damaging_missense:other_missense,pLoF:damaging_missense:other_missense:synonymous,pLoF:damaging_missense:other_missense:synonymous:non_coding
    # pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense,pLoF:damaging_missense:other_missense:synonymous
    # pLoF,damaging_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:synonymous

    dx-upload-all-outputs
}
