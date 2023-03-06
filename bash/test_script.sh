#!/usr/bin/env bash
#
# @description Write Hail Table of variants 
# @author Frederik Heymann Lassen (with edits by Nik Baya)
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=write_variants_ht
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas
#SBATCH --output=logs/write_variants_ht.log
#SBATCH --error=logs/write_variants_ht.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --requeue

sleep 100
