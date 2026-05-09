#!/usr/bin/env bash

GENOTYPES="$1" # Exome sequencing PLINK bed file (*.bed) or bgen (*.bgen) prefix
ANC="$2" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
PHENO_GROUP="$3"
PRED="$4" # .list file created by REGENIE step 1
PHENOCOL="$5" # Phenotype column
PHENOFILE_DNAX=$6
COVARFILE_DNAX=$7
flags="$8"
MAF_OR_AAF="$9"
ANNOT_TEMPLATE=${10} # where <FILETYPE> in the template is replaced with "setlist" or "annotations"
CHROM="${11}" # Chromosome
OUT="${12}" # Output prefix


# List of individuals who satisfy both:
# - Pass WES QC and are in ancestry group
# - Pass genotype array QC
KEEP="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes.qced.${ANC}.id"


# Create updated pred .list file
PRED_LOCAL="$HOME/tmp-predfile.txt"
cat ${PRED} | sed 's/\/home\/dnanexus\/out\/out\///g' > ${PRED_LOCAL} # Remove output prefix from .loco file
PRED=${PRED_LOCAL}
head $PRED

PHENOFILE_LOCAL="$HOME/tmp-phenofile.tsv"
if [ "$PHENO_GROUP" = "Height" ] || [[ "$PHENO_GROUP" == "original_phenos"* ]] || [[ "$PHENO_GROUP" == "microalbumin_urine_qced" ]] || [[ "${PHENO_GROUP}" == "standard_prs_controls" ]] || [[ "${PHENO_GROUP}" == "standardprs_"*"covariateresid_v2_2"* ]] || [[ "${PHENO_GROUP}" == "prs_as_covariate_v2_2_"* ]]; then
  gunzip -c $PHENOFILE_DNAX > $PHENOFILE_LOCAL
elif [[ "$PHENO_GROUP" == "mahalanobis_v2_"* ]]; then
  # Add FID field
  gunzip -c $PHENOFILE_DNAX | awk '{ print $1,$0 }' | sed '1 s/^IID/FID/g' > $PHENOFILE_LOCAL
else
  echo "ERROR: Unrecognised PHENO_GROUP $PHENO_GROUP. Could not process PHENOFILE_DNAX into PHENOFILE_LOCAL." && exit 1
fi

echo "PHENOFILE_LOCAL"
head $PHENOFILE_LOCAL

COVARFILE_LOCAL="$HOME/tmp-covarfile.tsv"
if [[ "${PHENOFILE_DNAX}" != "${COVARFILE_DNAX}" ]]; then
  gunzip -c $COVARFILE_DNAX > $COVARFILE_LOCAL
else
  ln -s ${PHENOFILE_LOCAL} ${COVARFILE_LOCAL}
fi

ANNOT_DIR="$( dirname ${ANNOT_TEMPLATE} )"
readonly ANNO="$( echo $ANNOT_TEMPLATE | sed 's/FILETYPE/annotations/g' )"
readonly SETLIST="$( echo $ANNOT_TEMPLATE | sed 's/FILETYPE/setlist/g' )"
readonly MASK="${ANNOT_DIR}/regenie_masks.txt"

# Define genotypes flag
if [ ${GENOTYPES} == *.bed ]; then
  BFILE=$( echo $GENOTYPES | sed 's/.bed$//g' )

  # Rename FID column in PLINK bfile
  awk '{ print $2,$2,$3,$4,$5,$6 }' ${BFILE}.fam > ${BFILE}.fam-tmp
  mv ${BFILE}.fam-tmp ${BFILE}.fam
  head ${BFILE}.fam
  
  # Edit .bim files to avoid REGENIE error when sets include multiple chroms
  if [[ "${CHROM}" == "GROUPED" ]]; then
    # Use dummy chromosome "23"
    awk '{ print 23,$2,$3,$4,$5,$6 }' ${BFILE}.bim > ${BFILE}.bim-tmp
    mv ${BFILE}.bim-tmp ${BFILE}.bim
    head ${BFILE}.bim
  fi

  genotypes_flag="--bed ${BFILE}"
elif [ ${GENOTYPES} == *.bgen ]; then
  BGEN=${GENOTYPES}
  SAMPLE=$( echo $GENOTYPES | sed 's/.bgen$/.sample/g' )
  genotypes_flag="--bgen ${BGEN} --sample ${SAMPLE}"
fi

get_allele_freq_file() {
  # PLINK .frq file format:
  # Column 2: 'SNP': Variant identifier
  # Column 5: 'MAF': Allele 1 frequency
  # NOTE: What PLINK reports as A1 in the .frq file is not the A1 in the bim file. Instead, it is the minor allele, such that the MAF column in .frq is indeed the minor allele frequency, not the allele frequency of allele 1 in the bim file.
  # NOTE: grep -v "NA" excludes variants where no individuals had defined genotypes and thus MAF couldn't be calculated
  _chrom=$1
  echo "/mnt/project/nbaya/regenie/data/allele_freq/ukb_wes_450k.qced.chr${_chrom}.${ANC}.${PHENO_GROUP}.${PHENOCOL}.frq"
}

get_command_to_read_uncompressed_or_compressed_file() {
  _filename=$1

  if [ -f "${_filename}" ]; then
    echo "cat ${_filename}"
  elif [ -f "${_filename}.gz" ]; then
    echo "gunzip -c ${_filename}"
  else
    echo "Error: No files matching $_filename{,.gz} were found" >&2
    exit 1
  fi
}


if [[ "$MAF_OR_AAF" == "MAF" ]]; then
  # AAF file to be constructed from PLINK .frq files
  AAF_FILE="${HOME}/tmp-aaf_file.txt"

  if [[ "$CHROM" =~ ^([1-9]|1[0-9]|2[0-2]|X)$ ]]; then 
    # If chrom in {1..22} or X
    # Get MAF from PLINK output and use as 'AAF' in place of actual AAF
    allele_freq="$( get_allele_freq_file ${CHROM} )"
    input_command="$( get_command_to_read_uncompressed_or_compressed_file $allele_freq )"
    eval ${input_command} | awk '{ print $2,$5 }' | tail -n+2 | grep -v "NA" > ${AAF_FILE}
  else
    # If not single chrom (e.g. grouped genes)
    # Extract chroms from ANNO file (no "chr" prefix, e.g. "12" instead of "chr12")
    # Assumes that 1st column of ANNO file contains variant IDs starting with "chr<CHROM>:"
    chroms=($( cat ${ANNO} | awk -F':' '{ print $1 }' | sed 's/^chr//g' | sort | uniq ))
    
    # Combine frq files for extracted chroms into a single file
    for chrom in ${chroms[@]}; do 
      tmp_allele_freq="$( get_allele_freq_file ${chrom} )"
      input_command=$( get_command_to_read_uncompressed_or_compressed_file $tmp_allele_freq )
      eval ${input_command} | awk '{ print $2,$5 }' | tail -n+2 | grep -v "NA" >> ${AAF_FILE}
    done
  fi

  aaf_file_flag="--aaf-file ${AAF_FILE}"
elif [[ "$MAF_OR_AAF" == "AAF" ]]; then
  aaf_file_flag="" # No need to provide file, allele freqs are calculated on the fly
else
  echo "ERROR: Unrecognised MAF_OR_AAF=${MAF_OR_AAF}. Could not define aaf_file_flag." && exit 1
fi


regenie \
  --step 2 \
  ${genotypes_flag} \
  --phenoFile ${PHENOFILE_LOCAL} \
  --covarFile ${COVARFILE_LOCAL} \
  ${flags} \
  --keep ${KEEP} \
  --pred ${PRED} \
  --anno-file ${ANNO} \
  --set-list ${SETLIST} \
  --mask-def ${MASK} \
  --vc-tests "skat,skato,acato" \
  ${aaf_file_flag} \
  --minMAC 0.5 \
  --bsize 400 \
  --maxCatLevels 25 \
  --out $OUT
  
if ! compgen -G "*regenie" > /dev/null; then
  echo "Error: No files matching *regenie found."
  exit 1
fi

mv *regenie ${OUT}.regenie

gzip *regenie

# Include -f to avoid error if no "tmp-" files present
rm -f tmp-*
