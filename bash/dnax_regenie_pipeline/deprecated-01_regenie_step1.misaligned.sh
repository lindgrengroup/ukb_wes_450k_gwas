#!/usr/bin/env bash

ANC="$1" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
GWAS_ID="$2" #"mahalanobis_v2_1_irnt_alpha0.01_upper_vs_inlier"
#PHENOCOL="$2" # Phenotype column name (e.g. "Height")


#ORIGINAL_BFILE="ukb22418_b0_v2.autosomes"

#BED_DNAX="nbaya/regenie/data/genotypes/$ORIGINAL_BFILE"
#BED_LOCAL=$ORIGINAL_BFILE
#dx download $BED_DNAX.{bed,bim,fam} --overwrite --output $HOME/
#BFILE_LOCAL=$3 # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")
BFILE_LOCAL="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes"


#PHENOFILE_DNAX="/mnt/project/nbaya/regenie/data/phenotypes/ukb.standing_height.20250508.tsv.gz"
PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.$GWAS_ID.tsv.gz"
PHENOFILE_LOCAL="$HOME/tmp-phenofile.tsv"
#gunzip -c $PHENOFILE_DNAX > $PHENOFILE_LOCAL
gunzip -c $PHENOFILE_DNAX | awk '{ print $1,$0 }' | sed '1 s/^IID/FID/g' > $PHENOFILE_LOCAL
head $PHENOFILE_LOCAL

# Convert phenotype list into comma-delimited list
# Exclude dichotomous case-control phenotypes by filtering out rows with "pearsonresid"
PHENOCOLLIST=$( cat /mnt/project/saige_pipeline/data/phenotypes/phenotype_list.$GWAS_ID.txt | grep -v "pearsonresid" | paste -sd ',' )
echo $PHENOCOLLIST

COVARFILE_DNAX="/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
COVARFILE_LOCAL="$HOME/tmp-covarfile.tsv"
gunzip -c $COVARFILE_DNAX > $COVARFILE_LOCAL

#COVARCOLLIST="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
#CATEGCOVARCOLLIST="sex"
COVARCOLLIST="sequencing_tranche" # Only covariate to include because all others have been regressed out
CATEGCOVARCOLLIST="sequencing_tranche"

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
# --bt : Trait is binary (0=control, 1=case, NA=missing)
# --ref-first : Reference allele is first in the PLINK file
trait_flag="--bt --ref-first" 
OUT="regenie_step1_${ANC}_${GWAS_ID}" #PHENOCOL}"

# Parameters based on default suggested by Regenie 'UKBB Analysis' recommendations:
# https://rgcgithub.github.io/regenie/recommendations/#step-1
#  --phenoCol ""$PHENOCOL"" \
regenie \
  --step 1 \
  --bed "$BFILE_LOCAL" \
  --phenoFile "$PHENOFILE_LOCAL" \
  --phenoColList ""${PHENOCOLLIST}"" \
  --covarFile "$COVARFILE_LOCAL" \
  --covarColList ""${COVARCOLLIST}"" \
  --catCovarList=""${CATEGCOVARCOLLIST}"" \
  --keep "$KEEP" \
  --extract "$EXTRACT" \
  --threads=${n_threads} \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp_lowmem \
  ${trait_flag} \
  --out $OUT

#rm tmp-*
