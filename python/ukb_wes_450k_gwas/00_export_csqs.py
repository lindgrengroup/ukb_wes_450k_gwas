#!/usr/bin/env python3
#
# @author Frederik Heymann Lassen (with edits by Nik Baya)
# Updated: 2022-11-04

import hail as hl
import argparse

from ukb_utils import hail_init
from ko_utils import csqs_case_builder


def main(args):

    in_vcf = args.in_vcf
    by = args.by
    by_explode = args.by_explode
    out_prefix = args.out_prefix
    out_type = args.out_type
    skip_and_write_ht_of_variants = bool(args.skip_and_write_ht_of_variants)

    hail_init.hail_bmrc_init(
        log='logs/hail/get_csqs.log',
        default_reference='GRCh38'
    )
    #hl._set_flags(no_whole_stage_codegen='1')


    mt = hl.import_vcf(
        path=in_vcf,
        force_bgz=True,
        reference_genome='GRCh38'
    )

    if skip_and_write_ht_of_variants:
        rows = mt.rows().select()
        rows.write(f'{out_prefix}.ht', overwrite=True)
        return

    if by_explode:
        mt = mt.explode_rows(mt.consequence.vep[by])
    mt = mt.annotate_rows(
        consequence_category=csqs_case_builder(
            worst_csq_expr=mt.consequence.vep[by],
            use_loftee=True
        )
    )
    ht = mt.rows()
    varid = hl.delimit([
        hl.str(ht.locus.contig),
        hl.str(ht.locus.position),
        ht.alleles[0],
        ht.alleles[1]],
        ':'
    )
    csqs = ht.consequence.vep[by]
    ht = ht.annotate(
        varid=varid,
        csqs=csqs
    )
    ht = ht.select('rsid', 'info', 'MAF', 'MAC', 'varid',
                   'csqs', 'consequence_category')
    # export as hail table or flattened tsv
    if out_type == "ht":
        ht.write(out_prefix + ".ht")
    else:
        ht.flatten().export(out_prefix + ".tsv.gz")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--in_vcf', default=None,
                        required=True, help='Path to input')
    parser.add_argument('--by', default=None, required=False,
                        help='What should be used for gene annotation')
    parser.add_argument('--by_explode', action='store_true',
                        help='What should be used for gene annotation')
    parser.add_argument('--out_prefix', default=None,
                        required=True, help='Path prefix for output dataset')
    parser.add_argument('--out_type', default=None,
                        required=False, help='either "ht" or "tsv".')
    parser.add_argument('--skip_and_write_ht_of_variants', default=False,
                        action='store_true', help='Whether to skip consequence annotation and just write a Hail table of the variants (locus + alleles)')

    args = parser.parse_args()

    main(args)
