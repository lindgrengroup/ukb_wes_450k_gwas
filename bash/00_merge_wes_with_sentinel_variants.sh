#!/bin/bash

set -e # Stop job if any command fails

## [OPTION] chrom
## Which chromosome to perform merge on
chrom=$1


## [OPTION] all_sentinel_variants
## File containing all sentinel variants, across all phenotypes, across all chromosomes
## NOTE: Rerunning finemapping was completed on 2023-09-18 for loci with windows that were incorrectly merged. So any versions following that date also have that change compared to previous versions.
#all_sentinel_variants="/mnt/project/saige_pipeline/data/sentinel_variants/sentinel.log10bf_gt_2.obesity.tsv.gz" # DEPRECATED: This list of sentinels was defined as variants with log10bf > 2
all_sentinel_variants="/mnt/project/saige_pipeline/data/sentinel_variants/sentinels.max_log10bf.include_ties.obesity.tsv.gz" # DEFAULT (2023-09-19): Sentinels defined as variants with maximum log10bf in a finemapped locus. We include all variants which are tied for maximum log10bf in given locus.

# Get sentinel variants for chromosome
extract_fname="tmp-sentinel_variants.rsids.chr${chrom}.txt"
gunzip -c ${all_sentinel_variants} | awk -v chrom=$chrom -v pheno=$pheno '{ if ($3==chrom) print $2 }' > ${extract_fname}
echo "Unique rsids to be extracted: $( uniq ${extract_fname} | wc -l )"

original_bgen_prefix="/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chrom}_b0_v3"
sentinel_bgen="tmp-sentinel_variants.chr${chrom}.bgen"

# Filter bgen to sentinel variants
bgenix \
    -g "${original_bgen_prefix}.bgen" \
    -incl-rsids ${extract_fname} \
    > ${sentinel_bgen}

# Create bgen index
bgenix -index -g ${sentinel_bgen}

# Convert from bgen to PLINK
sentinel_bfile="tmp-sentinel_variants.chr${chrom}"

plink2 \
--bgen "${sentinel_bgen}" ref-first \
--sample "${original_bgen_prefix}.sample" \
--make-bed \
--out ${sentinel_bfile}

wes="/mnt/project/Barney/wes/sample_filtered/ukb_wes_450k.qced.chr${chrom}"

# Set FIDs of sentinel bfile to 0 to match WES bfile
awk '{ print 0,$2,$3,$4,$5,$6 }' ${sentinel_bfile}.fam | column -t > tmp-${sentinel_bfile}.fam
mv tmp-${sentinel_bfile}.fam ${sentinel_bfile}.fam

# Merge
merged="ukb_wes_450k.qced.merged_w_imputed_sentinels.chr${chrom}"

plink \
  --bfile ${wes} \
  --bmerge ${sentinel_bfile} \
  --make-bed \
  --out ${merged}

# Clean up files (all remaining files are uploaded as output)
rm ${merged}.nosex
rm tmp-*

#for suffix in {bed,bim,fam}; do
#    dx upload ./${merged}.${suffix} --path "project-GBvkP10Jg8Jpy18FPjPByv29:/saige_pipeline/data/sentinel_variants/wes_merged_w_sentinels/"
#done
#
#rm ./*


