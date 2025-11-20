#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np


finemapping_dir = '/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_finemapping'
annot_dir = '/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/annotations/imputed_v3'


def get_ensembl_transcripts():
    ensembl_transcripts = pd.read_csv(
        '/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas/data/annotations/imputed_v3/ensembl_transcripts.GRCh37.tsv.gz',
        sep='\t',
        compression='gzip'
    )
    
    return ensembl_transcripts


def get_ensembl_genes():
    ensembl_transcripts = get_ensembl_transcripts()
    # NOTE: We only care about genes, not transcripts (of which several may compose a gene)
    # Therefore, we will only keep data relevant at the gene level
    
    # NOTE: Ensembl uses 1-based coordinate system:
    # https://www.ensembl.org/Help/Faq?id=286#:~:text=Ensembl%20uses%20a%20one%2Dbased,UCSC%20use%20the%20same%20genome.
    genes = ensembl_transcripts[['ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position']].drop_duplicates()

    return genes
    

def get_ensembl_gene_id_to_symbol_dict(ensembl_genes=None):
    if ensembl_genes is None:
        genes = get_ensembl_genes()
    gene_id_to_symbol_dict = dict(zip(*genes[['ensembl_gene_id','hgnc_symbol']].T.values))

    return gene_id_to_symbol_dict


def get_imputed_v3_variants(chrom):
    # NOTE: It seems that the imputed v3 variants exported from the original bgen files
    # have the same indexing as Ensembl (and therefore must be 1-indexed):
    # http://grch37.ensembl.org/Homo_sapiens/Variation/Explore?db=core;r=22:16050749-16051749;v=rs62224609;vdb=variation;vf=383210739
    variants = pd.read_csv(
        f'{finemapping_dir}/association_processing/bgen_v1.1.4-CentOS6.8-x86_64/imputed_v3.pass_qc.variant_ids.chr{chrom}.txt',
        sep='\t',
        skiprows=1, # skip comment at top of file
        skipfooter=1, # skip comment at the bottom of file
        dtype={
            'chromosome_name': str,
        }
    )
    variants = variants.drop(columns={'alternate_ids'})
    assert all(variants['number_of_alleles']==2), 'Not all variants are biallelic'
    return variants


def map_snp_to_gene(variant_position):
    genes_overlapping_variant = genes_in_chrom[
        (genes_in_chrom['start_position'] <= variant_position)
        & (variant_position <= genes_in_chrom['end_position'])
    ]['ensembl_gene_id']
    
    return ';'.join(genes_overlapping_variant)


def get_is_overlapping(gene_start, gene_stop, variant_start, variant_stop):
    overlap_start = max(gene_start, variant_start) 
    overlap_stop = min(gene_stop, variant_stop)
    return overlap_start <= overlap_stop


def map_indel_to_genes(variant_start, variant_stop):
    is_overlapping = genes_in_chrom.apply(
        lambda row: get_is_overlapping(
            gene_start=row['start_position'], 
            gene_stop=row['end_position'],
            variant_start=variant_start,
            variant_stop=variant_stop
        ),
        axis=1
    )
    genes_overlapping_variant = genes_in_chrom[is_overlapping]['ensembl_gene_id']
    
    return ';'.join(genes_overlapping_variant)


def get_imputed_v3_gene_annotations(chrom):
    variants = get_imputed_v3_variants(chrom=chrom)
    ensembl_genes = get_ensembl_genes()
    genes_in_chrom = ensembl_genes[ensembl_genes['chromosome_name']==str(chrom)]

    gene_id_to_symbol_dict = get_ensembl_gene_id_to_symbol_dict(ensembl_genes=ensembl_genes)

    # Annotate SNPs
    is_snp = (
        (variants['first_allele'].str.len()==1)
        & (variants['alternative_alleles'].str.len()==1)
    )
    variants.loc[is_snp,'ensembl_gene_id'] = variants.loc[is_snp,'position'].apply(map_snp_to_gene)

    # Annotate indels
    is_indel = ~is_snp
    variants.loc[is_indel,'ensembl_gene_id'] = variants
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', help='chromosome to run')
    args = parser.parse_args()

    get_imputed_v3_gene_annotations(chrom=args.chrom)
