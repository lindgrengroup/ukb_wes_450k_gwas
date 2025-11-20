#!/usr/bin/env python3

import argparse
import pandas as pd

from ukb_wes_450k_gwas.utils import get_obesity_phenotype_list

OBESITY_PHENOS = get_obesity_phenotype_list().tolist()

UKB_WES_450K_GWAS_DIR = '/gpfs3/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas'


def get_imputed_v3_covs():
    covs = pd.read_csv(
        f'{UKB_WES_450K_GWAS_DIR}/data/covariates/covariates.tsv.gz',
        sep='\t',
        compression='gzip'
    )
    return covs


def get_wes_450k_obesity_phenos_and_covs(sex="both_sexes"):
    phenos = pd.read_csv(
        f'{UKB_WES_450K_GWAS_DIR}/data/phenotypes/ukb_obesity_phenos_and_covariates.{sex}.tsv.gz',
        compression='gzip',
        sep='\t'
    )
    return phenos


def get_imputed_v3_obesity_phenos_and_covs_path(sex="both_sexes"):
    return f'{UKB_WES_450K_GWAS_DIR}/data/phenotypes/ukb_imputed_v3_obesity_phenos_and_covariates.{sex}.tsv'


def write_imputed_v3_obesity_phenos_and_covs(sex="both_sexes"):
    covs = get_imputed_v3_covs()

    phenos = get_wes_450k_obesity_phenos_and_covs(sex=sex)
    phenos = phenos[['FID', 'IID']+OBESITY_PHENOS]

    merge = covs.merge(phenos, on=['FID', 'IID'])

    # Note: SAIGE requires file to be uncompressed
    path = get_imputed_v3_obesity_phenos_and_covs_path(sex=sex)
    merge.to_csv(
        path,
        sep='\t',
        index=False
    )

    print(f'File written to {path}')


def get_imputed_v3_obesity_phenos_and_covs():
    return pd.from_csv(
        get_imputed_v3_obesity_phenos_and_covs(),
        sep='\t'
    )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--sex', default="both_sexes",
                        help='Sex to use. Options: "both_sexes", "female", "male"')
    args = parser.parse_args()
    sex = args.sex
    write_imputed_v3_obesity_phenos_and_covs(sex=sex)
