#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="02_regenie_group_test.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}
#dx rm -f ${script_dnax}
#dx upload ${script_local} --path ${script_dnax}

# OPTIONS
readonly PHENOCOL="Height" # Phenotype column name (e.g. "Height")
anc="EUR" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
#for anc in {"AFR","EAS"}; do
  
  PRED="regenie_step1_${anc}_${PHENOCOL}_pred.list"
  LOCO="regenie_step1_${anc}_${PHENOCOL}_1.loco"
  STEP1_DIR="nbaya/regenie/data/step1/${anc}"
  PRED_PATH="/mnt/project/${STEP1_DIR}/${PRED}"
  LOCO_PATH="${STEP1_DIR}/${LOCO}"
  
  for chrom in {14,23,20}; do
    if [ $chrom -eq 23 ]; then
      chrom="X"
    fi
  
    bfile="ukb_wes_450k.qced.chr${chrom}"
    bfile_path="/Barney/wes/sample_filtered/${bfile}"
  
    echo "Running REGENIE group test with ${bfile} for $anc $PHENOCOL"
  
    #if [ $chrom -eq 1 ] || [ $chrom -eq 19 ]; then
    #  instance_type="mem1_ssd2_v2_x8"
    #else
    #  instance_type="mem1_ssd1_v2_x8"
    #fi
  	#plink_bfile="/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr${chrom}"
  
    bed_size=$( dx ls -l "${bfile_path}.bed" 2> /dev/null | cut -f5 -d' ' )
    if (( $( echo "${bed_size} > 240" | bc ) )); then
      instance_type="mem1_ssd1_v2_x36"
      priority="high"
    elif (( $( echo "${bed_size} > 120" | bc ) )); then
      instance_type="mem1_ssd1_v2_x16"
      priority="high"
    elif (( $( echo "${bed_size} > 50" | bc ) )); then
      instance_type="mem1_ssd1_v2_x8"
      #instance_type="mem2_ssd1_v2_x8" && echo "############ OVERRIDING TEMPORARILY: mem3_ssd1_v2_x4 machines on low priority are fickle, using mem2_ssd1_v2_x8 instead "
      #priority="low"
      priority="high"
    else
      instance_type="mem1_ssd1_v2_x4"
      priority="high"
    fi
    
    out="regenie_group_test.${anc}.${PHENOCOL}.chr${chrom}"
  
    destination="nbaya/regenie/data/step2/group_tests/${anc}"
    dx mkdir -p ${destination}
  
    dx run swiss-army-knife \
    	-iin="${script_dnax}" \
      -iin="${bfile_path}.bed" \
      -iin="${bfile_path}.bim" \
      -iin="${bfile_path}.fam" \
      -iin=${LOCO_PATH} \
    	-icmd="bash ${script} ${bfile}.bed ${anc} ${PRED_PATH} ${PHENOCOL} ${chrom} ${out}" \
    	--name="regenie_group_test_${anc}_${PHENOCOL}_c${chrom}" \
    	--instance-type "$instance_type" \
    	--priority="$priority" \
    	--destination="$destination" \
    	--brief \
    	-y
  
  done
#done
