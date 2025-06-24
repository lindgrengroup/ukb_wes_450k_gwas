#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="02_regenie_variant_test.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}
#dx rm -f ${script_dnax}
#dx upload ${script_local} --path ${script_dnax}

# OPTIONS
readonly PHENOCOL="Height" # Phenotype column name (e.g. "Height")

#readonly anc="EUR" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
for anc in {"EAS","SAS"}; do
  anc_lower=$( echo ${anc} | tr '[:upper:]' '[:lower:]' )
  
  PRED="regenie_step1_${anc}_${PHENOCOL}_pred.list"
  LOCO="regenie_step1_${anc}_${PHENOCOL}_1.loco"
  STEP1_DIR="/nbaya/regenie/data/step1/${anc}"
  PRED_PATH="/mnt/project/${STEP1_DIR}/${PRED}"
  LOCO_PATH="${STEP1_DIR}/${LOCO}"
  
  # OPTION 1: Pseudovariants (additive, recessive encodings for various consequence masks)
  #BFILE_DIR="/wes_ko_ukbb/data/phased/encode_alt_qced/${anc_lower}/spliceai=0.50_cadd=28.1_revel=0.773/vcf_plus_plink"
  ##csq="pLoF_damaging_missense" # Options: damaging_missense_lc, other_missense, nonsynonymous, synonymous, damaging_missense, pLoF_damaging_missense, pLoF
  #for csq in {synonymous,nonsynonymous,other_missense,damaging_missense,pLoF}; do
  ##encoding="additive" # Options: additive, recessive
  #for encoding in {recessive,additive}; do
  #bfile="UKB.wes.merged.phased.full_qc.${anc_lower}.af05.popmax.pp0.90.spliceai=0.50_cadd=28.1_revel=0.773.${csq}.${encoding}.auto"
  #bfile_path="${BFILE_DIR}/${bfile}"
  ##bgen="UKB.wes.merged.phased.full_qc.eur.af05.popmax.pp0.90.spliceai=0.50_cadd=28.1_revel=0.773.pLoF_damaging_missense.recessive.auto.bgen"
  ##bgen_path="/nbaya/regenie/data/genotypes/${bgen}"
  ##sample_path=$( echo $bgen_path | sed 's/.bgen$/.sample/g' )
  #out="regenie_variant_test.${anc}.${PHENOCOL}.${csq}.${encoding}"


  # OPTION 2: All variants
  for chrom in {1..23}; do
    if [ $chrom -eq 23 ]; then
      chrom="X"
    fi

  bfile="ukb_wes_450k.qced.chr${chrom}"
  bfile_path="/Barney/wes/sample_filtered/${bfile}"
  out="regenie_variant_test.${anc}.${PHENOCOL}.chr${chrom}"
  
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
  
  destination="/nbaya/regenie/data/step2/variant_tests/${anc}"
  dx mkdir -p ${destination}
  
  echo "Running REGENIE step 2 (out=${out})"
  
  dx run swiss-army-knife \
  	-iin="${script_dnax}" \
    -iin="${bfile_path}.bed" \
    -iin="${bfile_path}.bim" \
    -iin="${bfile_path}.fam" \
    -iin="${LOCO_PATH}" \
  	-icmd="bash ${script} ${bfile}.bed ${anc} ${PRED_PATH} ${PHENOCOL} ${out}" \
  	--name="${out}" \
  	--instance-type "$instance_type" \
  	--priority="$priority" \
  	--destination="${destination}" \
  	--brief \
  	-y

  #break

done
done
#done
