#!/usr/bin/env bash

GENOTYPES="$1" # Exome sequencing PLINK bed file (*.bed) or bgen (*.bgen) prefix
ANC="$2" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
pheno_group="$3"
PRED="$4" # .list file created by REGENIE step 1
PHENOCOL="$5" # Phenotype column
PHENOFILE_DNAX=$6
COVARFILE_DNAX=$7
flags="$8"
MAF_OR_AAF="$9"
ANNOT_VERSION=${10} # Either "v6" or "v7"
CHROM="${11}" # Chromosome
OUT="${12}" # Output prefix

# # --- START DEBUG BLOCK ---
# echo "--- Checking arguments ---"
# echo "7 (COVARFILE_DNAX): $7"
# echo "8 (flags): $8"
# echo "9 (MAF_OR_AAF): $9"
# echo "10 (ANNOT_VERSION): ${10}"
# echo "11 (CHROM): ${11}"
# echo "12 (OUT): ${12}"
# echo "--------------------------"
# exit 1 # Stop the script from running

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
if [ "$pheno_group" = "Height" ] || [[ "$pheno_group" == "original_phenos"* ]] || [[ "$pheno_group" == "microalbumin_urine_qced" ]] || [[ "${pheno_group}" == "standard_prs_controls" ]] || [[ "${pheno_group}" == "standardprs_"*"covariateresid_v2_2"* ]]; then
  gunzip -c $PHENOFILE_DNAX > $PHENOFILE_LOCAL
elif [[ "$pheno_group" == "mahalanobis_v2_"* ]]; then
  # Add FID field
  gunzip -c $PHENOFILE_DNAX | awk '{ print $1,$0 }' | sed '1 s/^IID/FID/g' > $PHENOFILE_LOCAL
else
  echo "ERROR: Unrecognised pheno_group $pheno_group. Could not process PHENOFILE_DNAX into PHENOFILE_LOCAL." && exit 1
fi
echo "PHENOFILE_LOCAL"
head $PHENOFILE_LOCAL

COVARFILE_LOCAL="$HOME/tmp-covarfile.tsv"
if [[ "${PHENOFILE_DNAX}" != "${COVARFILE_DNAX}" ]]; then
  gunzip -c $COVARFILE_DNAX > $COVARFILE_LOCAL
else
  ln -s ${PHENOFILE_LOCAL} ${COVARFILE_LOCAL}
fi

ANNOT_DIR="/mnt/project/nbaya/regenie/data/annotations/${ANNOT_VERSION}"
readonly ANNO="${ANNOT_DIR}/regenie_annotations.${ANNOT_VERSION}.chr${CHROM}.txt"
readonly SETLIST="${ANNOT_DIR}/regenie_setlist.${ANNOT_VERSION}.chr${CHROM}.txt"
readonly MASK="${ANNOT_DIR}/regenie_masks.txt"

# Define genotypes flag
if [ ${GENOTYPES} == *.bed ]; then
  BFILE=$( echo $GENOTYPES | sed 's/.bed$//g' )
	
	# Rename FID column in PLINK bfile
  awk '{ print $2,$2,$3,$4,$5,$6 }' ${BFILE}.fam > ${BFILE}.fam-tmp
  mv ${BFILE}.fam-tmp ${BFILE}.fam
  head ${BFILE}.fam

  genotypes_flag="--bed ${BFILE}"
elif [ ${GENOTYPES} == *.bgen ]; then
  BGEN=${GENOTYPES}
  SAMPLE=$( echo $GENOTYPES | sed 's/.bgen$/.sample/g' )
  genotypes_flag="--bgen ${BGEN} --sample ${SAMPLE}"
fi

if [[ "$MAF_OR_AAF" == "MAF" ]]; then
  # Get MAF from PLINK output and use as 'AAF' in place of actual AAF
  allele_freq="/mnt/project/nbaya/regenie/data/allele_freq/ukb_wes_450k.qced.chr${CHROM}.${ANC}.${pheno_group}.${PHENOCOL}.frq"

  # Read uncompressed or compressed version
  if [ -f "${allele_freq}" ]; then
    input_command="cat ${allele_freq}"
  elif [ -f "${allele_freq}.gz" ]; then
    input_command="gunzip -c ${allele_freq}"
  else
    echo "ERROR: No files matching ${allele_freq}{,.gz} were found"
    exit 1
  fi
  AAF_FILE="${HOME}/tmp-aaf_file.txt"

  # PLINK .frq file format:
  # Column 2: 'SNP': Variant identifier
  # Column 5: 'MAF': Allele 1 frequency
  # NOTE: What PLINK reports as A1 in the .frq file is not the A1 in the bim file. Instead, it is the minor allele, such that the MAF column in .frq is indeed the minor allele frequency, not the allele frequency of allele 1 in the bim file.
  # NOTE: grep -v "NA" excludes variants where no individuals had defined genotypes and thus MAF couldn't be calculated
  eval ${input_command} | awk '{ print $2,$5 }' | tail -n+2 | grep -v "NA" > ${AAF_FILE}
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

mv *regenie ${OUT}.regenie

gzip *regenie

# Include -f to avoid error if no "tmp-" files present
rm -f tmp-*
