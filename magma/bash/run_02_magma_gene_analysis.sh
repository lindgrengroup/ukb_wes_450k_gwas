#!/usr/bin/env bash
#
# Run MAGMA gene-level analysis
#
# Author: Nik Baya (2023-06-05)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=02_magma
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/magma_gene_analysis.log
#SBATCH --error=logs/magma_gene_analysis.errors.log
#SBATCH --open-mode=append
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-23
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
chrom=$( get_chr ${SLURM_ARRAY_TASK_ID} )

pheno_list=( $( cat "${GWAS_WD}/data/phenotypes/phenotype_list.${phenotype_group}.txt" ) )
pheno="${pheno_list[${pheno_idx}]}"
gwas_id="${pheno}-${pop}-${sex}"

# GWAS sumstats (pre-computed with SAIGE)
sumstats="${MAGMA_WD}/sumstats/saige_variant_test.${gwas_id}-imputed_v3.chr${chrom}.tsv"
n_samples=$(( $( head -2 ${sumstats} | tail -1 | awk '{ print $5 }' ) ))

# Gene locations (columns reordered to use HGNC symbols)
gene_loc="${MAGMA_WD}/gene_locations/NCBI37.3.gene.loc.reformat"

# Gene annotation results (chromosome specific)
gene_annot="${MAGMA_WD}/output/${phenotype_group}/${pop}/${sex}/${gwas_id}.chr${chrom}.genes.annot"

# Reference genotypes
# bfile="${MAGMA_WD}/ref_data/g1000_eur" # 1000 Genomes Europeans
bfile="${GWAS_WD}/depict/data/genotype_data_plink/ukbb.fat-distn.bolt.association.imputed.hrc.geno_0.01.info_0.8.depict.unrelateds"
#bfile="${GWAS_WD}/data/00_set_up/ukb_impv3.subset_to_wes_450k_qc_pass_eur.chr${chrom}"

## [OUTPUT]
out_dir="${MAGMA_WD}/output/${phenotype_group}/${pop}/${sex}"
out="${out_dir}/${gwas_id}-ukbrefpanel"
mkdir -p "${out_dir}"


${magma} \
  --bfile ${bfile} \
  --pval ${sumstats} N=${n_samples} \
  --gene-annot ${gene_annot} \
  --batch ${chrom} "chr" \
  --out ${out}

