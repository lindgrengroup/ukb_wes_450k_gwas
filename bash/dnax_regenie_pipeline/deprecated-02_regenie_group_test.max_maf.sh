#!/usr/bin/env bash

GENOTYPES="$1" # Exome sequencing PLINK bed file (*.bed) or bgen (*.bgen) prefix
ANC=$2 # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
PRED="$3" # .list file created by REGENIE step 1
PHENOCOL="$4" # Phenotype column
MAX_MAF="$5" # Maximum minor allele frequency for inclusion of variants in burden mask
CHROM="$6" # Chromosome
OUT=$7 # Output prefix

# List of individuals who satisfy both:
# - Pass WES QC and are in ancestry group
# - Pass genotype array QC
KEEP="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes.qced.${ANC}.id"

# Remove output prefix from .loco file
PRED_LOCAL="$HOME/tmp-predfile.txt"
cat ${PRED} | sed 's/\/home\/dnanexus\/out\/out\///g' > ${PRED_LOCAL}
PRED=${PRED_LOCAL}
head $PRED

GWAS_ID=$PHENOCOL
if [ "$GWAS_ID" = "Height" ]; then
  PHENOFILE_DNAX="/mnt/project/nbaya/regenie/data/phenotypes/ukb.standing_height.20250508.tsv.gz"
  trait_flag="--qt --apply-rint --ref-first"
elif [ "$GWAS_ID" = "mahalanobis_v2_1_irnt_alpha0.01"* ]; then
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.$GWAS_ID.tsv.gz"
  trait_flag="--bt --ref-first"
elif [ "$GWAS_ID" = "original_phenos_qt" ]; then
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.original_phenos.tsv.gz"
  trait_flag="--qt --apply-rint --ref-first"
elif [ "$GWAS_ID" = "microalbumin_urine_qced" ]; then
    PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.qced_biomarkers.tsv.gz"
    trait_flag="--qt --apply-rint --ref-first"
else
  echo "ERROR: Unrecognised GWAS_ID $GWAS_ID. Could not define PHENOFILE_DNAX." && exit 1
fi
PHENOFILE_LOCAL="$HOME/tmp-phenofile.tsv"
gunzip -c $PHENOFILE_DNAX > $PHENOFILE_LOCAL
head $PHENOFILE_LOCAL

COVARFILE_DNAX="/mnt/project/nbaya/regenie/data/phenotypes/ukb_brava_default_covariates.20250508.tsv.gz"
COVARFILE_LOCAL="$HOME/tmp-covarfile.tsv"
gunzip -c $COVARFILE_DNAX > $COVARFILE_LOCAL

ANNOT_DIR="/mnt/project/nbaya/regenie/data/annotations/v7"
readonly ANNO="${ANNOT_DIR}/regenie_annotations.chr${CHROM}.txt"
readonly SETLIST="${ANNOT_DIR}/regenie_setlist.chr${CHROM}.txt"
readonly MASK="${ANNOT_DIR}/regenie_masks.txt"

# NOTE: This was never run
#ANNOT_VERSION=6
#ANNOT_DIR="/mnt/project/nbaya/regenie/data/annotations/v${ANNOT_VERSION}"
#readonly ANNO="${ANNOT_DIR}/regenie_annotations.v${ANNOT_VERSION}.chr${CHROM}.txt"
#readonly SETLIST="${ANNOT_DIR}/regenie_setlist.v${ANNOT_VERSION}.chr${CHROM}.txt"
#readonly MASK="${ANNOT_DIR}/regenie_masks.pLoF_damaging_missense.txt" # Only includes masks for pLoF-only and damaging missense only

# Get MAF from PLINK output and use as 'AAF' in place of actual AAF

allele_freq="/mnt/project/nbaya/regenie/data/allele_freq/ukb_wes_450k.qced.chr${CHROM}.${ANC}.${PHENOCOL}.frq"
AAF_FILE="${HOME}/tmp-aaf_file.txt"
awk '{ print $2,$5 }' ${allele_freq} | tail -n+2 | grep -v "NA" > ${AAF_FILE}


trait_flag="--qt --apply-rint --phenoCol=${PHENOCOL}"
#trait_flag="--bt --firth --approx --pThresh 0.1"

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


regenie \
  --step 2 \
  ${genotypes_flag} \
  --phenoFile ${PHENOFILE_LOCAL} \
  --covarFile ${COVARFILE_LOCAL} \
  ${trait_flag} \
  --keep ${KEEP} \
  --pred ${PRED} \
  --anno-file ${ANNO} \
  --set-list ${SETLIST} \
  --mask-def ${MASK} \
  --aaf-file ${AAF_FILE} \
  --aaf-bins ${MAX_MAF} \
  --vc-maxAAF ${MAX_MAF} \
  --vc-tests "skat,skato,acato" \
  --minMAC 0.5 \
  --bsize 400 \
  --out $OUT

mv *regenie ${OUT}.regenie

gzip *regenie

#rm tmp-*
