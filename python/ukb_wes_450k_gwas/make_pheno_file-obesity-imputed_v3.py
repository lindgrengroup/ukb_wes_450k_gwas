#!/usr/bin/env python3

import pandas as pd

UKB_WES_450K_GWAS_DIR='/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas'

def get_imputed_v3_covs():
    covs = pd.read_csv(
            f'{UKB_WES_450K_GWAS_DIR}/data/covariates/covariates.tsv.gz', 
            sep='\t', 
            compression='gzip'
    )
    return covs

def get_wes_450k_obesity_phenos_and_covs():
    phenos = pd.read_csv(
            f'{UKB_WES_450K_GWAS_DIR}/data/phenotypes/ukb_obesity_phenos_and_covariates.both_sexes.tsv', 
            sep='\t'
    )
    return phenos

def get_imputed_v3_obesity_phenos_and_covs_path():
    return f'{UKB_WES_450K_GWAS_DIR}/data/phenotypes/ukb_imputed_v3_obesity_phenos_and_covariates.both_sexes.tsv'

def write_imputed_v3_obesity_phenos_and_covs():
    covs = get_imputed_v3_covs()

    phenos = get_wes_450k_obesity_phenos_and_covs()
    phenos = phenos[[c for c in phenos.columns if c not in phenos.columns[2:30]]]
    assert len(phenos.columns)==12 # Check that columns are FID, IID + 10 obesity phenotypes

    merge = covs.merge(phenos, on=['FID', 'IID'])

    merge.to_csv(
            get_imputed_v3_obesity_phenos_and_covs_path(),
            sep='\t',
            index=False
            )

    print(f'File written to {get_imputed_v3_obesity_phenos_and_covs_path()}')

def get_imputed_v3_obesity_phenos_and_covs():
    return pd.from_csv(
            get_imputed_v3_obesity_phenos_and_covs(),
            sep='\t'
            )

if __name__=='__main__':
    write_imputed_v3_obesity_phenos_and_covs()
