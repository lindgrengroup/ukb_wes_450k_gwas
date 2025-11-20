#!/bin/bash

set -u # throws error if variables are undefined

WD="/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/bash/dnax_regenie_pipeline"

readonly script="01_regenie_step1.sh"
#readonly script="01_regenie_step1.misaligned.sh"

script_local="${WD}/${script}"
script_dnax="/nbaya/regenie/bash/${script}"

source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file ${script_local} ${script_dnax}
#dx rm -f ${script_dnax}
#dx upload ${script_local} --path ${script_dnax}

## OPTIONS
# readonly pheno_group="Height" 
#readonly pheno_group="mahalanobis_v2_2_irnt_upper_vs_inlier" # Aligned status (upper vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals
#readonly pheno_group="mahalanobis_v2_2_irnt_lower_vs_inlier" # Aligned status (lower vs. inlier); v2.2 uses unrelated individuals; v2.1 uses all available individuals 
#readonly pheno_group="original_phenos_qt" # Quantitative traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="original_phenos_bt" # Dichotomous traits (the original phenotypes corresponding to those used in misaligned project)
#readonly pheno_group="standard_prs_controls" 
# WARNING:standardprs_covariateresid_v2_2 should be run in cases and controls separately - otherwise, the mean imputation in REGENIE step 1 can be messed up due to cases and controls being mutatually exclusive.
#readonly pheno_group="standardprs_covariateresid_v2_2_ctrls" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
#readonly pheno_group="standardprs_covariateresid_v2_2_cases" # NOT SUBSET TO UNRELATED, despite saying v2.2. Uses PRS z-score with all GWAS covariates regressed out, except for sequencing tranche.
readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_ctrls" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche
#readonly pheno_group="standardprs_unrelated_covariateresid_v2_2_cases" # Subset to unrelateds separately in cases and controls. Uses PRS (not z-score of PRS), then residualised for GWAS covariates except for sequencing tranche
#readonly pheno_group="microalbumin_urine_qced"

# Get variables corresponding to pheno_group
export PHENO_GROUP=$pheno_group
source _load_variables_for_pheno_group.sh # This script requires a global variable PHENO_GROUP


# Define efine phenotype file and flags relevant to phenotype typeancestry group
readonly anc="EUR" # Genetic ancestry group (options: "EUR", "AFR", "EAS", "SAS", "AMR")
#for anc in {"EAS","SAS","AMR"}; do
  #readonly BFILE="ukb22418_b0_v2.autosomes" # PLINK bfile to use to fit model (e.g. "ukb22418_b0_v2.autosomes")
  
  echo "Running REGENIE step 1 on $anc for ${pheno_group}"
  
  destination="/nbaya/regenie/data/step1/${anc}"
  dx mkdir -p ${destination}
  
  dx run swiss-army-knife \
  	-iin="${script_dnax}" \
  	-icmd="""bash ${script} $anc $pheno_group ${PHENOFILE_DNAX} \"${pheno_col_flag}\" $COVARFILE_DNAX \"${covar_flag}\" \"${trait_flag}\" """ \
  	--name="regenie_step1_${anc}_${pheno_group}" \
  	--instance-type "mem1_ssd1_v2_x8" \
  	--priority="high" \
  	--destination="${destination}" \
  	--brief \
  	-y

  
#done
