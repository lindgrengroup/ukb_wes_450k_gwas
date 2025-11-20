#!/usr/bin/env bash
#
# Run MAGMA gene set enrichment analysis
#
# Author: Nik Baya (2023-06-05)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=03_magma
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/magma_gene_set_analysis.log
#SBATCH --error=logs/magma_gene_set_analysis.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --requeue

#set -eu # Stop job if any command fails

readonly GWAS_WD="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas"
readonly MAGMA_WD="${GWAS_WD}/magma"
source "/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/resources/ukb_utils/bash/cluster_utils.sh"

## [EXECUTABLE]
magma="${MAGMA_WD}/magma_static/magma"

## [INPUT]
phenotype_group="obesity"
pop="eur"
pheno_idx=$1 # 0-indexed
sex=$2

pheno_list=( $( cat "${GWAS_WD}/data/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
pheno="${pheno_list[${pheno_idx}]}"
gwas_id="${pheno}-${pop}-${sex}"

# Set annotation
# Gene set to use
# Options: cp, go
gene_set="cp"

if [ "${gene_set}" = "go" ]; then
  # Gene Ontology (GO) gene sets with HGNC symbols (n=16k)
  set_annot="${MAGMA_WD}/gene_sets/c5.all.v2023.1.Hs.symbols.gmt"
elif [ "${gene_set}" = "cp" ]; then 
  # Canonical pathways gene sets with HGNC symbols (n=3090)
  set_annot="${MAGMA_WD}/gene_sets/c2.cp.v2023.1.Hs.symbols.gmt"
fi

## [OUTPUT]
out_prefix="${MAGMA_WD}/output/${phenotype_group}/${pop}/${sex}/${gwas_id}-ukbrefpanel"
gene_results="${out_prefix}.genes.raw"
out_set_analysis="${out_prefix}.geneset_${gene_set}"

# Create merged gene results file if it does not exist
if [ ! -f ${gene_results} ]; then
  ${magma} \
    --merge ${out_prefix} \
    --out ${out_prefix}
fi

${magma} \
  --gene-results ${gene_results} \
  --set-annot ${set_annot} \
  --out ${out_set_analysis}

