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
readonly PHENO_GROUP="Height" # Uses BRaVA default covariates.

# Get variables corresponding to `PHENO_GROUP`
source _load_variables_for_pheno_group.sh

# Build flags as an array. This is the most robust way to handle
# arguments that might contain spaces.
flags_array=()

# Add static flags. Assumes flags with values are passed as single items
# (e.g., trait_flag="--trait-file /path/to/trait")
[[ -n "${trait_flag}" ]] && flags_array+=(${trait_flag})
[[ -n "${covar_flag}" ]] && flags_array+=(${covar_flag})


#readonly anc="EAS" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
anc="EUR"
#for anc in {"EAS","SAS"}; do
  anc_lower=$( echo ${anc} | tr '[:upper:]' '[:lower:]' )
  
  STEP1_DIR="/nbaya/regenie/data/step1/${anc}"
  PRED="regenie_step1_${anc}_${PHENO_GROUP}_pred.list"
  PRED_PATH="${STEP1_DIR}/${PRED}"
  PRED_DNAX="/mnt/project/${STEP1_DIR}/${PRED}"

# pheno_idx: One-indexed (i.e. starts with 1, ends with n)
for pheno_idx in {1..1}; do

  # Get phenotype column name (first column in pred file, in which each row is a different phenotype)
  pheno_col=$( dx cat ${PRED_PATH} | sed "${pheno_idx}q;d" | awk '{ print $1 }')
  pheno_col_flag="--phenoCol=${pheno_col}"

  # Create a *copy* of the static flags array for this loop iteration
  local_flags_array=("${flags_array[@]}" "${pheno_col_flag}")

  # Convert array to a single string, joined by spaces.
  flags_string="${local_flags_array[*]}"

  # Define LOCO file
  LOCO="regenie_step1_${anc}_${pheno_col}_${pheno_idx}.loco"
  LOCO_PATH="${STEP1_DIR}/${LOCO}"

  # VARIANT_TYPE
  # - "pseudovariant_*" # Additive, recessive encodings for various consequence masks (for Fred/BRaVa)
  # - "regular"  # Normal variants, no encodings
  #VARIANT_TYPE="pseudovariant_original" # Initial version (Nik ran May 2025)
  VARIANT_TYPE="pseudovariant_chr5filtered" # Followup version just for chr5 plof+damaging_missense (Nik ran Dec 2025)
  #VARIANT_TYPE="regular"

  # Load job parameters into arrays
  bfile_path_array=()
  out_array=()

  if [[ "${VARIANT_TYPE}" == "pseudovariant_"* ]]; then

    # Encoding
    #encoding="additive" # Options: additive, recessive
    for encoding in {recessive,additive}; do

    # Consequence
    csq="pLoF_damaging_missense" # Options: damaging_missense_lc, other_missense, nonsynonymous, synonymous, damaging_missense, pLoF_damaging_missense, pLoF
    ##for csq in {synonymous,nonsynonymous,other_missense,damaging_missense,pLoF}; do

    if [[ "${VARIANT_TYPE}" == "pseudovariant_original" ]]; then
      bfile_dir="/wes_ko_ukbb/data/phased/encode_alt_qced/${anc_lower}/spliceai=0.50_cadd=28.1_revel=0.773/vcf_plus_plink"
      if [[ "${anc}" == "EUR" ]]; then
        bfile_dir+="/force_chr_name" # File path is changed for EUR only
      elif [[ "${anc}" == "AMR" ]]; then
        echo "Warning: No pseudovariants were created for anc=AMR. Skipping."
        continue
      fi
      bfile="UKB.wes.merged.phased.full_qc.${anc_lower}.af05.popmax.pp0.90.spliceai=0.50_cadd=28.1_revel=0.773.${csq}.${encoding}.auto"
      job_out="regenie_variant_test.${anc}.${pheno_col}.${csq}.${encoding}"
    elif [[ "${VARIANT_TYPE}" == "pseudovariant_chr5filtered" ]]; then
      if [[ "${anc}" != "EUR" ]] || [[ "${csq}" != "pLoF_damaging_missense" ]]; then
        echo "Warning: VARIANT_TYPE=pseudovariant_chr5filtered is only allowed for anc=EUR and csq=pLoF_damaging_missense. Skipping"
        continue
      fi
      bfile_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/eur/spliceai=0.50_cadd=28.1_revel=0.773/filtered_chr5/vcf_plus_plink/force_chr_name"
      bfile="UKB.wes.chr5.phased.qc_final.${anc_lower}.af05.popmax.pp0.90.spliceai=0.50_cadd=28.1_revel=0.773.${csq}.${encoding}.chr5.filtered"
      job_out="regenie_variant_test.${anc}.${pheno_col}.${csq}.${encoding}.chr5filtered"
    else
      echo "Error: Invalid VARIANT_TYPE ${VARIANT_TYPE}. Exiting."
      exit 1
    fi

    job_bfile_path="${bfile_dir}/${bfile}"
    bfile_path_array+=(${job_bfile_path})

    out_array+=(${job_out})

  done
  #done

  elif [[ "${VARIANT_TYPE}" == "regular" ]]; then
    
    for chrom in {1..23}; do
      if [ $chrom -eq 23 ]; then
        chrom="X"
      fi

      bfile="ukb_wes_450k.qced.chr${chrom}"
      job_bfile_path="/Barney/wes/sample_filtered/${bfile}"
      bfile_path_array+=(${job_bfile_path})

      job_out="regenie_variant_test.${anc}.${pheno_col}.chr${chrom}"
      out_array+=($job_out)
    done

  else
    echo "Error: Invalid VARIANT_TYPE ${VARIANT_TYPE}. Exiting."
    exit 1
  fi
  
  destination="/nbaya/regenie/data/step2/variant_tests/${PHENO_GROUP}/${anc}"
  dx mkdir -p ${destination}
  
  # Calculate number of jobs for this phenotype
  n_jobs=${#out_array[@]}

  for job_idx in `seq 0 $(( n_jobs-1 ))`; do

    # Load job parameters from arrays
    bfile_path=${bfile_path_array[$job_idx]}
    bfile="${bfile_path##*/}"
    out=${out_array[$job_idx]}

    # Get DNAnexus machine instance type based on the size of the bed file 
    instance_type=$( bash ./_get_instance_type_for_bed_file.sh "${bfile_path}.bed" ) 
    [[ $instance_type != "mem"* ]] && exit 1 # Exit if instance_type is invalid
    priority="high" # Set job priority (low/normal/high)

    echo "Running REGENIE step 2 (out=${out})"
    
    dx run swiss-army-knife \
    	-iin="${script_dnax}" \
      -iin="${bfile_path}.bed" \
      -iin="${bfile_path}.bim" \
      -iin="${bfile_path}.fam" \
      -iin="${LOCO_PATH}" \
    	-icmd="bash ${script} ${bfile}.bed ${anc} ${PHENO_GROUP} ${PRED_DNAX} ${pheno_col} ${PHENOFILE_DNAX} ${COVARFILE_DNAX} \"${flags_string}\" ${out}" \
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
#done
