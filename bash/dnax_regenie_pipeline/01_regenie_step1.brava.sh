#!/usr/bin/env bash

ANC="$1" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
PHENOCOL="$2" # Phenotype column name (e.g. "Height")


#ORIGINAL_BFILE="ukb22418_b0_v2.autosomes"

#BED_DNAX="nbaya/regenie/data/genotypes/$ORIGINAL_BFILE"
#BED_LOCAL=$ORIGINAL_BFILE
#dx download $BED_DNAX.{bed,bim,fam} --overwrite --output $HOME/
#BFILE_LOCAL=$3 # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")
BFILE_LOCAL="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes"


PHENOFILE_DNAX="/mnt/project/nbaya/regenie/data/phenotypes/ukb.standing_height.20250508.tsv.gz"
PHENOFILE_LOCAL="$HOME/tmp-phenofile.tsv"
gunzip -c $PHENOFILE_DNAX > $PHENOFILE_LOCAL

COVARFILE_DNAX="/mnt/project/nbaya/regenie/data/phenotypes/ukb_brava_default_covariates.20250508.tsv.gz"
COVARFILE_LOCAL="$HOME/tmp-covarfile.tsv"
gunzip -c $COVARFILE_DNAX > $COVARFILE_LOCAL

COVARCOLLIST="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
CATEGCOVARCOLLIST="sex"

# List of individuals who satisfy both:
# - Pass WES QC and are in ancestry group
# - Pass genotype array QC
KEEP="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes.qced.${ANC}.id"

#KEEP_WES_QC_DNAX="/mnt/project/brava/inputs/ancestry_sample_ids/qced_${anc}_sample_IDs.txt"
#KEEP_WES_QC_LOCAL="tmp-sample_ids.wes_qc_pass_$anc.txt"
#awk '{ print $1,$1 }' $KEEP_WES_QC_DNAX > $KEEP_WES_QC_LOCAL
#KEEP="$HOME/tmp-keep.txt"
#join <(sort ${KEEP_WES_QC_LOCAL}) <(sort ${KEEP_ARRAY_QC_DNAX}) | awk '{ print $1,$2 }' > $KEEP

# Define list of variants which passed genotype array QC
EXTRACT="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes.qced.${ANC}.snplist"
#EXTRACT_LOCAL="$HOME/tmp-extract.txt"
#cp $EXTRACT_DNAX $EXTRACT_LOCAL
#head -10000 $EXTRACT_DNAX > $EXTRACT_LOCAL

n_threads=16
trait_flag="--qt --apply-rint"
OUT="regenie_step1_${ANC}_${PHENOCOL}"

regenie \
  --step 1 \
  --bed "$BFILE_LOCAL" \
  --phenoFile "$PHENOFILE_LOCAL" \
  --phenoCol ""$PHENOCOL"" \
  --covarFile "$COVARFILE_LOCAL" \
  --covarColList ""${COVARCOLLIST}"" \
  --catCovarList=""${CATEGCOVARCOLLIST}"" \
  --keep "$KEEP" \
  --extract "$EXTRACT" \
  --threads=${n_threads} \
  --bsize 1000 \
  ${trait_flag} \
  --out $OUT

#rm tmp-*
