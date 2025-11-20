#!/usr/bin/env bash
#
# Run MAGMA SNP to gene mapping
#
# Author: Nik Baya (2023-06-05)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=01_magma
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/magma_annotate.log
#SBATCH --error=logs/magma_annotate.errors.log
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

# Gene locations (columns reordered to use HGNC symbols)
gene_loc="${MAGMA_WD}/gene_locations/NCBI37.3.gene.loc.reformat"

## [OUTPUT]
gwas_id="${pheno}-${pop}-${sex}"
out_dir="${MAGMA_WD}/output/${phenotype_group}/${pop}/${sex}"
out="${out_dir}/${gwas_id}.chr${chrom}"
mkdir -p ${out_dir}

## Reformat sumstats
original="${GWAS_WD}/data/02_saige_variant_test-imputed_v3/${phenotype_group}/${pop}/${sex}/saige_variant_test.${gwas_id}-imputed_v3.chr${chrom}.tsv.gz"
sumstats="${MAGMA_WD}/sumstats/saige_variant_test.${gwas_id}-imputed_v3.chr${chrom}.tsv"
#gunzip -c $original | awk '{ print $3,$1,$2,$13,$14 }' | head -1 > $sumstats
gunzip -c $original | awk '{ print $3,$1,$2,$13,$14 }'  >> $sumstats
sed -i 's/MarkerID/SNP/g; s/p.value/P/g' ${sumstats}


${magma} \
  --annotate \
  --snp-loc ${sumstats} \
  --gene-loc ${gene_loc} \
  --out ${out}
