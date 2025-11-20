#!/usr/bin/env python3

import argparse
import pandas as pd


def main(chrom):

    annots_dir = '/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/annotations/imputed_v3'

    # Get list of variant IDs from imputed v3 dataset
    variants = pd.read_csv(
        f'{annots_dir}/imputed_v3.pass_qc.variant_ids.chr{chrom}.txt', delim_whitespace=True, comment='#')
    variants['chrom:pos'] = 'chr' + \
        variants.chromosome.astype(str)+':'+variants.position.astype(str)

    
    annot = pd.read_csv(f'{annots_dir}/v2g_annots_hg19.chr{chrom}.txt', delim_whitespace=True,
                        names=['chrom', 'position', 'chrom:pos', 'gene_id:fpred_max_label', 'fpred_max_score'])

    annot = annot.drop_duplicates(subset=['chrom:pos','gene_id:fpred_max_label','fpred_max_score'])

    merged = variants.merge(annot[['chrom:pos', 'gene_id:fpred_max_label', 'fpred_max_score']], on=[
                            'chrom:pos'], how='inner')

    merged['ensembl_gene_id'] = merged['gene_id:fpred_max_label'].str.split(
        ':', expand=True)[0]
    merged['fpred_max_label'] = merged['gene_id:fpred_max_label'].str.split(
        ':', expand=True)[1]

    bridge = pd.read_csv('/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz',
                         compression='gzip', delim_whitespace=True)
    bridge = bridge.drop_duplicates(subset=['hgnc_symbol', 'ensembl_gene_id'])

    merged_w_hgnc = merged.merge(
        bridge[['hgnc_symbol', 'ensembl_gene_id']], on='ensembl_gene_id')

    merged_w_hgnc = merged_w_hgnc[['chromosome','position','rsid','ensembl_gene_id','hgnc_symbol','fpred_max_label','fpred_max_score']]

    merged_w_hgnc = merged_w_hgnc.drop_duplicates(subset=['rsid','ensembl_gene_id'])

    def concat_genes(x):
        return '/'.join(x)

    merged_w_hgnc = merged_w_hgnc[['rsid','hgnc_symbol']].sort_values(by='hgnc_symbol').groupby(['rsid']).agg(concat_genes)

    merged_w_hgnc = merged_w_hgnc.reset_index()

    merged_w_hgnc.to_csv(
        f'{annots_dir}/v2g_annots_hg19.hgnc_symbol.v2.chr{chrom}.tsv.gz',
        compression='gzip',
        index=False,
        sep='\t'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', help='chromosome to run')
    args = parser.parse_args()

    main(chrom=args.chrom)
