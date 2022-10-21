#!/bin/bash

saige_data_dir="wes_450k:/saige_pipeline/data"

instance_type="mem1_ssd1_v2_x4"
traitType=quantitative
invNormalize=FALSE
phenoCol="Body_mass_index_(BMI)_combined_geneticeur_50k_InvNorm"
covariatesList="age,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15,pc16,pc17,pc18,pc19,pc20,pc21,age2,is_female_age,is_female_age2"
qCovarColList="is_female,assessment_centre,genotyping_batch"
sampleIDCol=IID
pheno_file="${saige_data_dir}/phenotypes/ukb_obesity_phenos_and_covariates.tsv.gz"
PLINK_for_vr="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.for_vr.chr21"
sparseGRMfile="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
sparseGRM_sample_file="${saige_data_dir}/00_set_up/ukb_array.wes_450k_qc_pass_eur.pruned_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
jobname="bmi_eur_50k"

workflow_id="workflow-GJ3vjb0Jg8JXZP8J5zjz6fgy"


dx run ${workflow_id} \
      -istage-common.phenofile="${pheno_file}" \
      -istage-common.bedfile="${PLINK_for_vr}.bed" \
      -istage-common.bimfile="${PLINK_for_vr}.bim" \
      -istage-common.famfile="${PLINK_for_vr}.fam" \
      -istage-common.spGRMfile="${sparseGRMfile}" \
      -istage-common.spGRMSamplefile="${sparseGRM_sample_file}" \
      -istage-common.output_prefix="${outputPrefix}" \
      -istage-common.phenoCol="${phenoCol}" \
      -istage-common.traitType="${traitType}" \
      -istage-common.covariatesList="${covariatesList}" \
      -istage-common.qCovarColList="${qCovarColList}" \
      -istage-common.sampleIDCol="${sampleIDCol}" \
      -istage-common.invNormalize="${invNormalize}" \
      --folder="${saige_data_dir}/01_fit_null" \
      --yes \
      --name=${jobname} \
      --instance-type=${instance_type}