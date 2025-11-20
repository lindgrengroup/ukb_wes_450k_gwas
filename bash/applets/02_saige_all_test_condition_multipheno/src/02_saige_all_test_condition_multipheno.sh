#!/bin/bash
# 02_saige_all_test_condition_multipheno 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    # mkdir -p out/{output,log}_files
    mkdir -p out/{gene_assoc,marker_list,variant_assoc,log_file}

    # # dx_wd=$( project-GBvkP10Jg8Jpy18FPjPByv29:/ )
    saige_data_dir="project-GBvkP10Jg8Jpy18FPjPByv29://saige_pipeline/data"
    output_dir="${saige_data_dir}/02_saige_all_test_condition${catevr}/${phenotype_group}/${pop}/${sex}"
    
    # [INPUT] SAIGE null model file directory
    fit_null_output_prefix="${saige_data_dir}/01_fit_null${catevr}/${phenotype_group}/${pop}/${sex}"

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
    #else
    #  # Remove "chr" prefix from chromosome column in bim file
    #  sed -i 's/^chr//g' "${WD}/bfile.bim"
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
      gwas_id="${pheno_col}-${pop}-${sex}"

      # [INPUT] Sentinel variants
      if [ "${pheno_col}" == "bmi" ]; then
      # NOTE: BMI was called "body_mass_index_bmi" when running the common variant finemapping
        imp_pheno="body_mass_index_bmi"
      else
        imp_pheno=${pheno_col}
      fi

      # Get variant IDs of sentinel variants
      # NOTE: These are variants from a common-variant GWAS (using imputed data) which have been finemapped.
      # NOTE: We define a "sentinel" variant as the variant with the highest log10 bayes factor for being causal in a given locus.
      # NOTE: As of 2023-09-19, variants which are tied for highest log10bf are all included. Previously I broke ties (higher z-score wins ties).
      # NOTE: Output as comma-delimited list (e.g. rs101,rs122,rs203 )
      # NOTE: Trailing commas in the sentinel variant list do not seem to affect the results. Regardless, I've removed them for a cleaner look.
      all_sentinel_variants="${WD}/in/all_sentinel_variants/"*
      sentinel_list=$( gunzip -c ${all_sentinel_variants} | awk -v chrom=$chrom -v imp_pheno=$imp_pheno '{ ORS = "," }{ if ($3==chrom && $18==imp_pheno) print $2 }' | sed 's/,$//g' )

      # If there are no sentinel variants to condition on
      if [ "${sentinel_list}" == "" ]; then
        echo "No sentinels for ${imp_pheno} chr${chrom}"
        fake_log_file="no_sentinels.${pheno_col}.chr${chrom}.log"
        if [ $( dx ls -l ${output_dir}/${fake_log_file} 2> /dev/null | wc -l ) -eq 0 ]; then
          touch "out/log_file/${fake_log_file}"
        fi
        continue
      fi

      # [OUTPUT] Output files
      # NOTE (2023-09-22): Some of the jobs were failing due to excess memory requests (sometimes close to 2 Tb of RAM requested)
      # To fix this, I reduced the number of tests we run. In particular, I believe that not running tests for the union of all
      # variant consequences (i.e. using all coding variants) and instead using the "minimal" set of tests to get the results we
      # need (pLoF and pLoF+damaging_missense) helped reduce the memory required. The full results will have file names starting 
      # with "saige_all_test", whereas the adapted minimal versions start with "minimal_tests...-". The default prefix is now 
      # "minimal_tests_w_syn", because it runs the minimal set of tests as well as a test of synonymous variants, in case we want 
      # to check if our tests are well calibrated.
      original_output_prefix="saige_all_test_condition${catevr}.${gwas_id}.chr${chrom}"
      output_prefix="minimal_tests_w_syn-${original_output_prefix}"
      original_results_file="${output_dir}/${original_output_prefix}.tsv.gz"
      results_file="${output_dir}/${output_prefix}.tsv.gz"

      # Check that neither the original "full" results exist nor the minimal test results
      if [ $( dx ls -l ${original_results_file} 2> /dev/null | wc -l  ) -eq 0 ] && [ $( dx ls -l ${results_file} 2> /dev/null | wc -l ) -eq 0 ]; then
  
        echo "Starting $gwas_id"
        
        # # [INPUT] SAIGE null model files
        model_file="${fit_null_output_prefix}/${gwas_id}.rda"
        variance_ratios="${fit_null_output_prefix}/${gwas_id}.varianceRatio.txt"
        dx download --overwrite "${model_file}" -o "${WD}/model_file.rda"
        dx download --overwrite "${variance_ratios}" -o "${WD}/variance_ratios.txt"


        # ## Run script
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
            --GMMATmodelFile ${HOME}/model_file.rda  \
            --varianceRatioFile ${HOME}/variance_ratios.txt  \
            --sparseGRMFile ${HOME}/in/sparse_grm/*  \
            --sparseGRMSampleIDFile ${HOME}/in/sparse_grm_samples/* \
            --groupFile ${HOME}/in/group_file/* \
            --annotation_in_groupTest "pLoF,pLoF:damaging_missense,synonymous" \
            --maxMAF_in_groupTest 0.001,0.01 \
            --is_output_markerList_in_groupTest TRUE \
            --is_fastTest TRUE \
            --condition "${sentinel_list}" \
            --SAIGEOutputFile ${HOME}/${output_prefix}.tsv 2>&1 | tee ${HOME}/${output_prefix}.log

      
        gzip "${output_prefix}.tsv" 
        mv "${output_prefix}.tsv.gz" "out/gene_assoc/" 
        gzip "${output_prefix}.tsv.markerList.txt" 
        mv "${output_prefix}.tsv.markerList.txt.gz" "out/marker_list/${output_prefix}.markerlist.txt.gz"
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

    #DEFAULT (all tests)
    #        --annotation_in_groupTest "pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous" \
    #        --maxMAF_in_groupTest 0.0001,0.001,0.01 \

    # --groupFile ${HOME}/group_file.txt \
    # pLoF,damaging_missense,synonymous,non_coding,pLoF:damaging_missense,pLoF:damaging_missense:synonymous,pLoF:damaging_missense:synonymous:non_coding
    # pLoF,damaging_missense,other_missense,synonymous,non_coding,pLoF:damaging_missense,pLoF:damaging_missense:other_missense,pLoF:damaging_missense:other_missense:synonymous,pLoF:damaging_missense:other_missense:synonymous:non_coding
    # pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense,pLoF:damaging_missense:other_missense:synonymous
    # pLoF,damaging_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:synonymous

    dx-upload-all-outputs
}
