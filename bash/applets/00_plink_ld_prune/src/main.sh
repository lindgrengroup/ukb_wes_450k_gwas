#!/bin/bash
#
# This app performs LD pruning on UKB genotype array data
# Before pruning, samples are subset to Europeans who pass WES 450k QC 
# Pruning is parallelized across chromosomes

main() {

  #
  # Debugging setup
  # ---------------
  # The -e flag causes bash to exit at any point if there is any error,
  # the -o pipefail flag tells bash to throw an error if it encounters an error within a pipeline,
  # the -x flag causes bash to output each line as it is executed
  #

  set -e -x -o pipefail

  #
  # SECTION: LD prune
  # --------------------------------------------------------
  #
  prune_jobs=()
  for chrom in {22..23}; do
    if [ "${chrom}" -eq 23 ]; then
      chrom="X"
    fi
    prune_jobs+=($(dx-jobutil-new-job -igenotype_array_bfile="${genotype_array_bfile}" -isamples_w_superpop="${samples_w_superpop}" -ichrom="${chrom}" ld_prune))
  done
  
  input_flags=()
  for job in "${prune_jobs[@]}"; do
    input_flags+=("-ipruned_chrom=${job}:pruned_chrom")
    input_flags+=("-ipruned_fam=${job}:pruned_fam")
    input_flags+=("-ipruned_bim=${job}:pruned_bim")
    input_flags+=("-ipruned_bed=${job}:pruned_bed")
  done

  merged_bfile="merged_pruned"
  merge_job=$(dx-jobutil-new-job "${input_flags[@]}" -imerged_bfile="${merged_bfile}" merge_bfiles)

  # head "${merge_job}:all_bims_name"

  dx-jobutil-add-output merged_ld_pruned "${merge_job}:all_bims" --class=array:jobref

}


# ####################################################################
# # This function will LD prune per-chromosome files
# #
# # Arguments:
# #   chrom: Specified chromosome for current job
# #   genotype_array_bfile
# #   samples_w_superpop
# ####################################################################
ld_prune() {

  echo "Chromosome to be pruned '${chrom}'"
  echo "PLINK bfile '${genotype_array_bfile}'"
  echo "samples_w_superpop '${samples_w_superpop}'"

  input_prefix="raw_c${chrom}"

  for suffix in {fam,bim,bed}; do
    dx download \
      "$( echo "${genotype_array_bfile}" | sed -e "s/#CHROM#/${chrom}/g" ).${suffix}" \
      --output "${input_prefix}.${suffix}"
  done

  dx download "${samples_w_superpop}" -o "samples_w_superpop.tsv"

  out_prefix="ld_pruned_c${chrom}"

  plink \
    --bfile "${input_prefix}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "samples_w_superpop.tsv" ) \
    --indep-pairwise 50 5 0.05 \
    --out "${out_prefix}"

  # Extract set of pruned variants and export to bfile
  plink \
    --bfile "${input_prefix}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "samples_w_superpop.tsv" ) \
    --extract "${out_prefix}.prune.in" \
    --make-bed \
    --out "${out_prefix}"

  dx-jobutil-add-output pruned_chrom "${chrom}" --class=string
  for suffix in {fam,bim,bed}; do
    out_file=$(dx upload "${input_prefix}.${suffix}" --brief)
    dx-jobutil-add-output "pruned_${suffix}" "${out_file}" --class=file
  done
}

####################################################################

####################################################################
merge_bfiles() {

  set -e -x -o pipefail

  echo "merged_bfile: ${merged_bfile}"

  merge_list="merge_list.txt"
  out="ukb_array.wes_450k_qc_pass_eur.pruned"

  for chrom in "${!pruned_chrom[@]}"; do
    dx download ${pruned_fam[$i]} -o "c${chrom}.fam"
    dx download ${pruned_bim[$i]} -o "c${chrom}.bim"
    dx download ${pruned_bed[$i]} -o "c${chrom}.bed"

    echo "c${chrom}" >> merge_list.txt
  done

  plink \
    --merge-list merge_list.txt \
    --make-bed \
    --out "${out}"
  
  #
  # Upload files
  # ------------
  #
  for suffix in {fam,bim,bed}; do
    dx-jobutil-add-output merged_ld_pruned "${out}.${suffix}" --class=array:file
  done
}
