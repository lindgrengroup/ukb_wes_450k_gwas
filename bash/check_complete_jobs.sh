#!/usr/bin/env bash

PHENO_LIST_FNAME="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/phenotypes/phenotype_list.obesity.txt"
WES_VARIANT_TEST_DIR="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/02_saige_variant_test"
IMPUTED_VARIANT_TEST_DIR="${WES_VARIANT_TEST_DIR}-imputed_v3"
WES_SET_TEST_DIR="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/02_saige_set_test"
WES_ALL_TEST_DIR="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/02_saige_all_test"

phenos=( $( cat ${PHENO_LIST_FNAME} ) )


echo ">> WES DATA <<"
for sex in {both_sexes,male,female}; do
  echo $sex
  for pheno_idx in {0..9}; do 
    pheno=${phenos[$pheno_idx]}
    ct=$( ls -1 ${WES_VARIANT_TEST_DIR}/obesity/eur/${sex}/saige_*test.${pheno}-*tsv.gz 2> /dev/null | wc -l )
    echo -e "* ${pheno}: ${ct}/23"
  done
done

echo
echo ">> IMPUTED DATA <<"
for sex in {both_sexes,male,female}; do
  echo $sex
  for pheno_idx in {0..9}; do 
    pheno=${phenos[$pheno_idx]}
    ct=$( ls -1 ${IMPUTED_VARIANT_TEST_DIR}/obesity/eur/${sex}/saige_*test.${pheno}-*tsv.gz 2> /dev/null | wc -l )
    echo -e "* ${pheno}: ${ct}/23"
  done
done

echo
echo ">> WES DATA (SET TEST) <<"
for sex in {both_sexes,male,female}; do
  echo $sex
  for pheno_idx in {0..9}; do 
    pheno=${phenos[$pheno_idx]}
    ct=$( ls -1 ${WES_SET_TEST_DIR}/obesity/eur/${sex}/saige_*test.${pheno}-*tsv.gz 2> /dev/null | wc -l )
    echo -e "* ${pheno}: ${ct}/23"
  done
done

echo ">> WES DATA (ALL TEST) <<"
for sex in {both_sexes,male,female}; do
  echo $sex
  for pheno_idx in {0..8}; do 
    pheno=${phenos[$pheno_idx]}
    variant_ct=$( ls -1 ${WES_ALL_TEST_DIR}/obesity/eur/${sex}/saige_*test.${pheno}-*singleAssoc.txt.gz 2> /dev/null | wc -l )
    gene_ct=$( ls -1 ${WES_ALL_TEST_DIR}/obesity/eur/${sex}/saige_*test.${pheno}-*.tsv.gz 2> /dev/null | grep -v "markerList" | wc -l )
    echo -e "* ${pheno}:"
    echo -e "  - gene:    ${gene_ct}/23"
    echo -e "  - variant: ${variant_ct}/23"  
  done
done
