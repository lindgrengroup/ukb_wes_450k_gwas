#!/usr/bin/env bash
# 
# Run ANNOVAR on imputed v3 variants for a given chromosome


WD="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/annotations/imputed_v3/annovar"
finemap_dir="/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_finemapping/"

chrom=$1

if [ $chrom -eq 23 ]; then
  chrom="X";
fi;     

# NOTE: Variant ID file is assumed to be output by bgenix -list, for example:
#   alternate_ids  rsid         chromosome  position    number_of_alleles  first_allele  alternative_alleles
#   .              rs537417111  21          9574185     2                  G             C
#   .              rs572757938  21          9660864     2                  G             A
#   .              rs562470929  21          9674879     2                  G             C
# WARNING: The following code only works if 'number_of_alleles' is 2 for all rows
variant_ids="${finemap_dir}/association_processing/bgen_v1.1.4-CentOS6.8-x86_64/imputed_v3.pass_qc.variant_ids.chr${chrom}.txt"
out="${WD}/imputed_v3.chr$chrom"
avinput="${WD}/imputed_v3.chr$chrom.avinput"

# 1) Create avinput file
#    Columns should be:
#    - Chromosome
#    - Variant start position (1-indexed)
#    - Variant end position (1-indexed)
#    - Reference allele
#    - Alternate allele
#    - [Optional] Additional variant-level info (e.g. rsid)
# 2) Run annotate_variation using the avinput file
#    - Use build hg19
#    - Use the default reference data in the /humandb directory provided by ANNOVAR
grep -v "#\|alternate" ${variant_ids} \
  | awk '{ print $3, $4, $4+length($6)-1, $6, $7, $2 }' \
  > ${avinput} \
&& ${WD}/annotate_variation.pl \
  -out ${out} \
  -build hg19 \
  --separate \
  ${avinput} \
  ${WD}/humandb/
