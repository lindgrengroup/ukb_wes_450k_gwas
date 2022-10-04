#!/bin/bash

instance_type="mem1_ssd1_v2_x4"
traitType=binary
invNormalize=FALSE
phenoCol=value
covariatesList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,age,age_sex,age2,age2_sex,sex
qCovarColList=sex
sampleIDCol=userId
pheno_file=phecode-250.2-both_sexes_wowithdrawl_with450kWES.tsv
PLINK_for_vr=ukb.EUR.for_grm.pruned.plink.forvr
sparseGRMfile=ukb.EUR_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRM_sample_file=ukb.EUR_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
jobname=phecode-250.2_step1

workflow_id=workflow-GB2gx6jJ6y3j4P2f2K60Pfxz


dx run ${workflow_id} \
      -istage-common.phenofile=saige:/SAIGE_GENE/phenotype/${pheno_file} \
      -istage-common.bedfile=saige:/SAIGE_GENE/genotype/${PLINK_for_vr}.bed \
      -istage-common.bimfile=saige:/SAIGE_GENE/genotype/${PLINK_for_vr}.bim \
      -istage-common.famfile=saige:/SAIGE_GENE/genotype/${PLINK_for_vr}.fam \
      -istage-common.spGRMfile=saige:/SAIGE_GENE/spGRM/${sparseGRMfile} \
      -istage-common.spGRMSamplefile=saige:/SAIGE_GENE/spGRM/${sparseGRM_sample_file} \
      -istage-common.output_prefix=${outputPrefix} \
      -istage-common.phenoCol=${phenoCol} \
      -istage-common.traitType=${traitType} \
      -istage-common.covariatesList=${covariatesList} \
      -istage-common.qCovarColList=${qCovarColList} \
      -istage-common.sampleIDCol=${sampleIDCol} \
      -istage-common.invNormalize=${invNormalize} \
      --folder saige:/SAIGE_GENE/step1_output/ \
      --yes \
      --name=${jobname} \
      --instance-type=${instance_type}