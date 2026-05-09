#!/bin/bash

echo "Loading variables for phenotype group '$PHENO_GROUP'..."

# Define phenotype file and flags relevant to phenotype type
# --bt : Trait is binary (0=control, 1=case, NA=missing)
# --qt : Trait is quantitative (continuous)
# Define phenotype column flag (either phenoCol or phenoColList)
if [ "$PHENO_GROUP" = "Height" ]; then
  PHENOFILE_DNAX="/mnt/project/nbaya/regenie/data/phenotypes/ukb.standing_height.20250508.tsv.gz"
  trait_flag="--qt --apply-rint"

	PHENOCOL="Height"
  pheno_col_flag="--phenoCol=${PHENOCOL}"

elif [[ "$PHENO_GROUP" == "mahalanobis_v2_"* ]]; then
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.$PHENO_GROUP.tsv.gz"
  trait_flag="--bt"

	# Convert phenotype list into comma-delimited list
  # Exclude dichotomous case-control phenotypes by filtering out rows with "pearsonresid"
  PHENOCOLLIST=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.$PHENO_GROUP.txt | grep -v "pearsonresid" | paste -sd ',' )
  pheno_col_flag="--phenoColList=${PHENOCOLLIST}"

elif [ "$PHENO_GROUP" = "original_phenos_qt" ]; then
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.original_phenos.tsv.gz"
  trait_flag="--qt --apply-rint"
	
	PHENOCOLLIST=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep -v "t2d\|cad\|osteoporosis" | paste -sd ',' )
  pheno_col_flag="--phenoColList=${PHENOCOLLIST}"

elif [ "$PHENO_GROUP" = "original_phenos_bt" ]; then
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.original_phenos.tsv.gz"
  trait_flag="--bt"
	
	PHENOCOLLIST=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.original_phenos.txt | grep "t2d\|cad\|osteoporosis" | paste -sd ',' )
  pheno_col_flag="--phenoColList=${PHENOCOLLIST}"

elif [[ "$PHENO_GROUP" == "prs_as_covariate_v2_2"* ]]; then
  # Get the phenotype which you want to regress out PRS from by taking the suffix of the PHENO_GROUP
  # e.g. prs_as_covariate_v2_2_qt_height -> height
  # e.g. prs_as_covariate_v2_2_bt_t2d -> t2d
  pheno_with_prs="${PHENO_GROUP##*_}"

  # Get the base phenotype group, with the phenotype suffix removed
  pheno_group_base="${PHENO_GROUP%_${pheno_with_prs}}" # either prs_as_covariate_v2_2_qt or prs_as_covariate_v2_2_bt

  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.${pheno_group_base}.tsv.gz"
  if [[ "${PHENO_GROUP}" == "prs_as_covariate_v2_2_qt"* ]]; then
    trait_flag="--qt --apply-rint"
  elif [[ "${PHENO_GROUP}" == "prs_as_covariate_v2_2_bt"* ]]; then
    trait_flag="--bt"
  else
    echo "ERROR: Expected ${PHENO_GROUP} to contain '_qt' or '_bt'" && exit 1
  fi
	
  pheno_col_flag="--phenoColList=${pheno_with_prs}"

elif [[ "${PHENO_GROUP}" == "standard_prs_controls" ]]; then
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.$PHENO_GROUP.tsv.gz"
  trait_flag="--qt" # IRNT already applied to the PRS

	PHENOCOLLIST=$( dx cat /saige_pipeline/data/phenotypes/phenotype_list.$PHENO_GROUP.txt | paste -sd ',' )
  pheno_col_flag="--phenoColList=${PHENOCOLLIST}"

elif  [[ ${PHENO_GROUP} == "standardprs_"*"covariateresid_v2_2_"* ]]; then
  # 1. Get the last part of the string (e.g., "ctrls" or "cases")
  type_suffix="${PHENO_GROUP##*_}"

  # 2. Get the base name (e.g., "standardprs_covariateresid_v2_2")
  file_base="${PHENO_GROUP%_${type_suffix}}"

	# 3. Get the grep pattern by removing the trailing 's' (e.g., "ctrl" or "case")
	grep_pattern="${type_suffix%s}"
	
	# 4. Construct the file path
	pheno_list_file="/saige_pipeline/data/phenotypes/phenotype_list.${file_base}.txt"

	# 5. Get variables
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.${file_base}.tsv.gz"
  trait_flag="--qt" # IRNT already applied to the PRS

	PHENOCOLLIST=$( dx cat "${pheno_list_file}" | grep "${grep_pattern}" | paste -sd ',' )
  pheno_col_flag="--phenoColList=${PHENOCOLLIST}"

elif [ "$PHENO_GROUP" = "microalbumin_urine_qced" ]; then
  PHENOFILE_DNAX="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.qced_biomarkers.tsv.gz"
  trait_flag="--qt --apply-rint"

	PHENOCOL="microalbumin_urine_qced"
  pheno_col_flag="--phenoCol=${PHENOCOL}"
else
  echo "ERROR: Unrecognised PHENO_GROUP $PHENO_GROUP." && exit 1
fi


# Define covariate file and covariate list
if [ "$PHENO_GROUP" = "Height" ]; then
  COVARFILE_DNAX="/mnt/project/nbaya/regenie/data/phenotypes/ukb_brava_default_covariates.20250508.tsv.gz"
  COVARCOLLIST="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
  CATEGCOVARCOLLIST="sex"
elif [[ "$PHENO_GROUP" == "mahalanobis_v2_"* ]] || [[ "${PHENO_GROUP}" == "standardprs"*"covariateresid_v2_2_"* ]]; then
  # NOTE: Unlike 'standard_prs_controls', standardprs_covariateresid_v2_2 already has all covariates regressed out (except for sequencing_tranche)
  COVARFILE_DNAX="/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
  COVARCOLLIST="sequencing_tranche" # Only covariate to include because all others have been regressed out
  CATEGCOVARCOLLIST="sequencing_tranche"
elif [[ "$PHENO_GROUP" == "original_phenos"* ]] || [[ "$PHENO_GROUP" == "microalbumin_urine_qced" ]] || [[ "${PHENO_GROUP}" == "standard_prs_controls" ]] || [[ "${PHENO_GROUP}" == "prs_as_covariate_v2_2_"* ]]; then
  #COVARFILE_DNAX="/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
  COVARFILE_DNAX=$PHENOFILE_DNAX
  COVARCOLLIST="age,age2,is_female,is_female_age,is_female_age2,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15,pc16,pc17,pc18,pc19,pc20,pc21,assessment_centre,sequencing_tranche"
  
  if [[ "${PHENO_GROUP}" == "prs_as_covariate_v2_2_"* ]]; then
    COVARCOLLIST+=",prs_${pheno_with_prs}" # Include the phenotype PRS as a covariate
  fi

  CATEGCOVARCOLLIST="is_female,assessment_centre,sequencing_tranche"
else
  echo "ERROR: Unrecognised PHENO_GROUP $PHENO_GROUP. Could not define COVARCOLLIST." && exit 1
fi
covar_flag="--covarColList=${COVARCOLLIST} --catCovarList=${CATEGCOVARCOLLIST}"
