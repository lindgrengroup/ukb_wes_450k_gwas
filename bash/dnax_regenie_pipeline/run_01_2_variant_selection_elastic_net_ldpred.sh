#!/bin/bash

# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------
# Mode: "elastic_net" or "ldpred2"
MODE="elastic_net"
#MODE="ldpred2"

# Set the target GENE here
GENE="$1"

# Random Seed: Set an integer (e.g. 42) or leave empty for random.
SEED=$2

#RADIUS=1000000
RADIUS=2000000

PHENO_FILE="/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.standardprs_unrelated_covariateresid_v2_2.tsv.gz"
GFF_PATH="/mnt/project/resources/ensembl/Homo_sapiens.GRCh37.87.gff3.gz"
DNAX_DIR="/mnt/project/nbaya/outliers/data/imputed_v3_condition_variants/condition_pgen"

# Determine PHENO_COL based on GENE mapping
case "$GENE" in
    "HNF1A"|"HNF4A"|"GCK")
        PHENO_COL="t2d_case_prs_resid"
        ;;
    "APOB"|"PCSK9"|"ANGPTL3")
        PHENO_COL="cad_ctrl_prs_resid"
        ;;
    "LDLR")
        PHENO_COL="cad_case_prs_resid"
        ;;
    "SLC30A8")
        PHENO_COL="t2d_ctrl_prs_resid"
        ;;
    "COL1A1"|"COL1A2")
        PHENO_COL="osteoporosis_case_prs_resid"
        ;;
    *)
        echo "Error: No phenotype mapping found for gene $GENE"
        exit 1
        ;;
esac

echo "Configuration: Gene=$GENE -> Phenotype=$PHENO_COL"


# Test Run: Set number of samples (e.g. 1000). Set to "" or 0 to disable (run full).
#TEST_SAMPLES=1000
TEST_SAMPLES="" # If null, do full run

# Output location
PROJECT_OUT_DIR="/nbaya/outliers/data/imputed_v3_condition_variants/lasso_elastic_net/"
# Location to store the script on DNAnexus
SCRIPT_REMOTE_DIR="/nbaya/outliers/scripts/"

# -------------------------------------------------------------------
# 1. UPLOAD R SCRIPT
# -------------------------------------------------------------------
R_SCRIPT_NAME="variant_selection_elastic_net_ldpred.R"
LOCAL_SCRIPT_PATH="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/R/${R_SCRIPT_NAME}"
REMOTE_SCRIPT_PATH="${SCRIPT_REMOTE_DIR}${R_SCRIPT_NAME}"

echo "Uploading R script to $REMOTE_SCRIPT_PATH ..."

# Load 'upload_file' function from dnax_utils.sh:
source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_qc/bash/dnax_utils.sh"

# usage: upload_file <local_path> <remote_path>
upload_file "$LOCAL_SCRIPT_PATH" "$REMOTE_SCRIPT_PATH"

# -------------------------------------------------------------------
# 2. CONSTRUCT COMMAND & JOB NAME
# -------------------------------------------------------------------

# Base Job Name
JOB_NAME="${MODE}_${GENE}_r${RADIUS}_${PHENO_COL}"

# Base R command
R_CMD="Rscript $R_SCRIPT_NAME \
        --gene $GENE \
        --radius $RADIUS \
        --mode $MODE \
        --pheno_file $PHENO_FILE \
        --pheno_col $PHENO_COL \
        --gff_path $GFF_PATH \
        --dnax_dir $DNAX_DIR \
        --out_dir ."

# Add test flag if configured
if [ -n "$TEST_SAMPLES" ] && [ "$TEST_SAMPLES" -gt 0 ]; then
    echo "Configuring TEST RUN with $TEST_SAMPLES samples."
    R_CMD="$R_CMD --test_n $TEST_SAMPLES"
    
    # Append Test info to Job Name
    JOB_NAME="${JOB_NAME}_TEST${TEST_SAMPLES}"
fi

# Add seed flag if configured
if [ -n "$SEED" ]; then
    echo "Configuring Random Seed: $SEED"
    R_CMD="$R_CMD --seed $SEED"
    
    # Append Seed info to Job Name
    JOB_NAME="${JOB_NAME}_seed${SEED}"
fi

# -------------------------------------------------------------------
# 3. SUBMIT JOB
# -------------------------------------------------------------------

echo "Submitting Job: $JOB_NAME"
echo "Destination: $PROJECT_OUT_DIR"

dx run swiss-army-knife \
   --name "$JOB_NAME" \
   -iin="$REMOTE_SCRIPT_PATH" \
   -icmd="
      # A. Install Dependencies
      sudo R -e 'install.packages(c(\"data.table\", \"optparse\", \"remotes\"), repos=\"https://cloud.r-project.org\")'
      sudo R -e 'install.packages(\"bigsnpr\", repos=\"https://cloud.r-project.org\")'

      # B. Run Analysis Script
      echo 'Running: $R_CMD'
      $R_CMD
   " \
   --instance-type "mem1_ssd1_v2_x16" \
   --destination "$PROJECT_OUT_DIR" \
   --priority "high" \
   --brief \
   --yes
