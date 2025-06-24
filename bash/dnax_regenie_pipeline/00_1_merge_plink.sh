#!/usr/bin/env bash

rm -f merge_list.autosomes.txt
for chrom in {1..22}; do 
  bfile="ukb22418_c${chrom}_b0_v2"
  ln -f -s "/mnt/project/Bulk/Genotype Results/Genotype calls/$bfile"* ./
  echo $bfile >> merge_list.autosomes.txt
done


plink2 \
  --pmerge-list merge_list.autosomes.txt bfile \
  --make-bed \
  --out ukb22418_b0_v2.autosomes
