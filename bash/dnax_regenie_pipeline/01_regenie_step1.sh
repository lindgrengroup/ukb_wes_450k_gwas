#!/usr/bin/env bash

anc="$1" # Genetic ancestry group (e.g. "EUR", "AFR", "EAS", etc.)
pheno_group="$2" #"mahalanobis_v2_1_irnt_alpha0.01_upper_vs_inlier"
PHENOFILE_DNAX="$3"
pheno_col_flag="$4"
COVARFILE_DNAX="$5"
covar_flag=$6
trait_flag="$7"

#ORIGINAL_BFILE="ukb22418_b0_v2.autosomes"

#BED_DNAX="nbaya/regenie/data/genotypes/$ORIGINAL_BFILE"
#BED_LOCAL=$ORIGINAL_BFILE
#dx download $BED_DNAX.{bed,bim,fam} --overwrite --output $HOME/
#BFILE_LOCAL=$3 # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")
# --ref-first : Reference allele is first in the genotype array PLINK file (I have checked that this is true for the first 9 genotype array variants, including a variant where the alt allele was not a minor allele, which suggests the variants aren't ordered by minor allele in the PLINK files)
BFILE_LOCAL="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes"


PHENOFILE_LOCAL="$HOME/tmp-phenofile.tsv"
if [ "$pheno_group" = "Height" ] || [[ "$pheno_group" == "original_phenos"* ]] || [[ "$pheno_group" == "microalbumin_urine_qced" ]] || [[ "${pheno_group}" == "standard_prs_controls" ]] || [[ "${pheno_group}" == "standardprs_"*"covariateresid_v2_2_"* ]]; then
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
gunzip -c $COVARFILE_DNAX > $COVARFILE_LOCAL
echo "COVARFILE_LOCAL"
head $COVARFILE_LOCAL



# List of individuals who satisfy both:
# - Pass WES QC and are in ancestry group
# - Pass genotype array QC
KEEP="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes.qced.${anc}.id"

#KEEP_WES_QC_DNAX="/mnt/project/brava/inputs/ancestry_sample_ids/qced_${anc}_sample_IDs.txt"
#KEEP_WES_QC_LOCAL="tmp-sample_ids.wes_qc_pass_$anc.txt"
#awk '{ print $1,$1 }' $KEEP_WES_QC_DNAX > $KEEP_WES_QC_LOCAL
#KEEP="$HOME/tmp-keep.txt"
#join <(sort ${KEEP_WES_QC_LOCAL}) <(sort ${KEEP_ARRAY_QC_DNAX}) | awk '{ print $1,$2 }' > $KEEP

# Define list of variants which passed genotype array QC
EXTRACT="/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes.qced.${anc}.snplist"
#EXTRACT_LOCAL="$HOME/tmp-extract.txt"
#cp $EXTRACT_DNAX $EXTRACT_LOCAL
#head -10000 $EXTRACT_DNAX > $EXTRACT_LOCAL

n_threads=16
OUT="regenie_step1_${anc}_${pheno_group}" #PHENOCOL}"

# Parameters based on default suggested by Regenie 'UKBB Analysis' recommendations:
# https://rgcgithub.github.io/regenie/recommendations/#step-1
# --ref-first : Reference allele is first in the genotype array PLINK file (I have checked that this is true for the first 9 genotype array variants, including a variant where the alt allele was not a minor allele, which suggests the variants aren't ordered by minor allele in the PLINK files)
# NOTE: --maxCatLevels needs to be manually increased because there are 22 assessment centres
regenie \
  --step 1 \
  --bed "$BFILE_LOCAL" \
  --phenoFile "$PHENOFILE_LOCAL" \
  "${pheno_col_flag}" \
  --ref-first \
  --covarFile "$COVARFILE_LOCAL" \
  ${covar_flag} \
  --maxCatLevels 25 \
  --keep "$KEEP" \
  --extract "$EXTRACT" \
  --threads=${n_threads} \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix tmp_lowmem \
  ${trait_flag} \
  --out $OUT

#rm tmp-*
