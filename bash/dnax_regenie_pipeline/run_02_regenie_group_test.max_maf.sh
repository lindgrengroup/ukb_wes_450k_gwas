#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="02_regenie_group_test.max_maf.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}
#dx rm -f ${script_dnax}
#dx upload ${script_local} --path ${script_dnax}

# OPTIONS
readonly PHENOCOL="Height" # Phenotype column name (e.g. "Height")
for max_maf in {0.01,0.001,0.0001}; do # Maximum allele frequency of variants included in mask
#max_maf=0.01
#for anc in {"EAS","SAS","AFR"}; do # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
anc="SAS"
#for anc in {"AFR","EAS"}; do

  
  PRED="regenie_step1_${anc}_${PHENOCOL}_pred.list"
  LOCO="regenie_step1_${anc}_${PHENOCOL}_1.loco"
  STEP1_DIR="nbaya/regenie/data/step1/${anc}"
  PRED_PATH="/mnt/project/${STEP1_DIR}/${PRED}"
  LOCO_PATH="${STEP1_DIR}/${LOCO}"
  
  for chrom in {1..23}; do
    if [ $chrom -eq 23 ]; then
      chrom="X"
    fi
  
    out="regenie_group_test.${anc}.${PHENOCOL}.maxmaf${max_maf}.chr${chrom}"
    destination="nbaya/regenie/data/step2/group_tests/${anc}"
    num_files=$(dx ls "${destination}/${out}.regenie.gz" | wc -l) # Check if output file already exists

    if [ "$num_files" -eq 0 ]; then

      echo "Running REGENIE group test (out=${out})"
      dx mkdir -p ${destination}

      bfile="ukb_wes_450k.qced.chr${chrom}"
      bfile_path="/Barney/wes/sample_filtered/${bfile}"
  
      bed_size=$( dx ls -l "${bfile_path}.bed" 2> /dev/null | cut -f5 -d' ' )
      if (( $( echo "${bed_size} > 240" | bc ) )); then
        instance_type="mem1_ssd1_v2_x36"
        priority="high"
      elif (( $( echo "${bed_size} > 120" | bc ) )); then
        instance_type="mem1_ssd1_v2_x16"
        instance_type="mem1_ssd1_v2_x36" # TEMPORARY OVERRIDE
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

      #priority="low"

      dx run swiss-army-knife \
      	-iin="${script_dnax}" \
        -iin="${bfile_path}.bed" \
        -iin="${bfile_path}.bim" \
        -iin="${bfile_path}.fam" \
        -iin=${LOCO_PATH} \
      	-icmd="bash ${script} ${bfile}.bed ${anc} ${PRED_PATH} ${PHENOCOL} ${max_maf} ${chrom} ${out}" \
      	--name="${out}" \
      	--instance-type "$instance_type" \
      	--priority="$priority" \
      	--destination="$destination" \
      	--brief \
      	-y
    else
      echo "Skipping REGENIE group test (file already exists: ${destination}/${out}.regenie.gz)"
    fi
  
  done
done
#done
