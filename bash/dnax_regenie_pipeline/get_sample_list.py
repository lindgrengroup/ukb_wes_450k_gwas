#!/usr/bin/env python3
import argparse
import pandas as pd


def main(gwas_id, pheno_col, anc_list, sex='both_sexes'):
    """
    Prepares IID lists for GWAS by filtering samples based on phenotype data,
    covariates, and ancestry-specific QC lists.

    Args:
        gwas_id (str): A unique identifier for the GWAS run.
        pheno_col (str): The name of the phenotype column to use.
        anc_list (str): A comma-separated string of ancestry groups.
    """
    chrom = 21  # Can be any chrom

    bfile = f"ukb_wes_450k.qced.chr{chrom}"
    bfile_path = f"/mnt/project/Barney/wes/sample_filtered/{bfile}"
    fam_path = f'{bfile_path}.fam'

    fam = pd.read_csv(fam_path, sep='\s+', names=[
                      'fid', 'IID', 'MID', 'PID', 'SEX', 'PHENO'])

    if pheno_col == 'Height':
        PHENOFILE_DNAX = "/mnt/project/nbaya/regenie/data/phenotypes/ukb.standing_height.20250508.tsv.gz"
        COVARFILE_DNAX = "/mnt/project/nbaya/regenie/data/phenotypes/ukb_brava_default_covariates.20250508.tsv.gz"
    elif pheno_col == 'microalbumin_urine_qced':
        PHENOFILE_DNAX = "/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.qced_biomarkers.tsv.gz"
        COVARFILE_DNAX = "/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
    elif gwas_id.startswith('original_phenos'):
        PHENOFILE_DNAX = f"/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.original_phenos.tsv.gz"
        COVARFILE_DNAX = "/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
    elif gwas_id.startswith('standardprs_covariateresid_v2_2'):
        PHENOFILE_DNAX = f"/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.standardprs_covariateresid_v2_2.tsv.gz"
        COVARFILE_DNAX = "/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
    elif gwas_id.startswith('prs_as_covariate_v2_2_qt'):
        PHENOFILE_DNAX = f"/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.prs_as_covariate_v2_2_qt.tsv.gz"
        COVARFILE_DNAX = "/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
    elif gwas_id.startswith('prs_as_covariate_v2_2_bt'):
        PHENOFILE_DNAX = f"/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.prs_as_covariate_v2_2_bt.tsv.gz"
        COVARFILE_DNAX = "/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"
    else:
        PHENOFILE_DNAX = f"/mnt/project/saige_pipeline/data/phenotypes/ukb_wes.phenos_and_covariates.{gwas_id}.tsv.gz"
        COVARFILE_DNAX = "/mnt/project/saige_pipeline/data/covariates/ukb_wes_standard_covs.tsv.gz"

    phenos = pd.read_csv(PHENOFILE_DNAX, sep='\t', dtype={'FID': str, 'IID': str})
    covs = pd.read_csv(COVARFILE_DNAX, sep='\t', dtype={'FID': str, 'IID': str})

    phenos = phenos[['IID', pheno_col]]

    merged = phenos.merge(covs, on=['IID'])

    if sex=='female':
        merged = merged[merged['is_female']==1]
    elif sex=='male':
        merged = merged[merged['is_female']==0]
    else:
        assert sex=='both_sexes' # No need to filter merged, just check that sex is 'both_sexes'

    for anc in anc_list.split(','):
        print(f"\nProcessing ancestry: {anc.upper()}") 
        KEEP=f"/mnt/project/nbaya/regenie/data/genotypes/ukb22418_b0_v2.autosomes.qced.{anc}.id"
        keep = pd.read_csv(KEEP, sep='\t', names=['FID','IID'], dtype={'FID':str, 'IID':str})
        is_all_notna = merged.notna().all(axis=1)
        is_pass_qc = merged['IID'].isin(keep['IID']) # Genotype array QC

        out = f'ukb_wes_450k.qced.{anc}.{gwas_id}.{pheno_col}.iids.txt'
        print(f'Number of samples: {(is_all_notna&is_pass_qc).sum()}')
        print(f'- Is all not NA:   {is_all_notna.sum()}')
        print(f'- Is pass QC:      {is_pass_qc.sum()}')
        out_df = merged.loc[is_all_notna & is_pass_qc, 'IID']
        out_df.to_csv(out, index=False, header=False)

if __name__ == '__main__':
    # Set up the command-line argument parser
    parser = argparse.ArgumentParser(description="Generate sample IID lists for GWAS analysis.")

    # Define the command-line arguments
    parser.add_argument('gwas_id', type=str, 
                        help='A unique identifier for the GWAS run (e.g., "UKB_WES_v1").')
    parser.add_argument('pheno_col', type=str, 
                        help='The name of the phenotype column in the phenotype file (e.g., "Height").')
    parser.add_argument('anc_list', type=str, 
                        help='A comma-separated string of ancestry groups to process (e.g., "EUR,AFR").')
    parser.add_argument('sex', type=str, default='both_sexes',
            help='Options: both_sexes, female, male')

    # Parse the arguments provided by the user
    args = parser.parse_args()

    # Call the main function with the parsed arguments
    main(args.gwas_id, args.pheno_col, args.anc_list, args.sex)
