import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from os.path import isfile
from datetime import datetime
from os.path import getmtime


pd.options.mode.chained_assignment = None  # default='warn'

# +
OBESITY_PHENOS = [
    'body_mass_index_bmi',
    'whr_adj_bmi',
    'bmi_impedance',
    'body_fat_percentage',
    'visceral_adipose_tissue_volume_vat',
    'abdominal_fat_ratio',
    'gynoid_tissue_fatp',
    'android_tissue_fatp',
    'total_tissue_fatp',
    'tissuefatp_androidgynoidratio'
]

# The four phenotypes with high sample size (>300k samples)
high_sample_size_obesity_phenos = [
    'body_mass_index_bmi',
    'whr_adj_bmi',
    'bmi_impedance',
    'body_fat_percentage'
]

obesity_phenotype_to_label_dict = {
    'bmi': 'BMI',
    'body_mass_index_bmi': 'BMI',
    'whr_adj_bmi': 'WHRadjBMI',
    'bmi_impedance': 'Impedance BMI',
    'body_fat_percentage': 'Body fat percentage',
    'visceral_adipose_tissue_volume_vat': 'Visceral adipose tissue',
    'abdominal_fat_ratio': 'Abdominal fat ratio',
    'gynoid_tissue_fatp': 'Gynoid tissue fat percentage',
    'android_tissue_fatp': 'Android tissue fat percentage',
    'total_tissue_fatp': 'Total tissue fat percentage',
    'tissuefatp_androidgynoidratio': 'Android-gynoid fat percentage ratio'
}
# -


SET_TEST_GROUPS = [
    'pLoF',
    'damaging_missense',
    'synonymous',
    'pLoF;damaging_missense',
    'pLoF;damaging_missense;synonymous',
    'Cauchy'
]

GWS_PVAL_THRESHOLD = 5e-8

CURRENT_CONSEQUENCE_VERSION = '6'

# +
# pheno_to_label_dict = dict(zip(BIOMARKERS, BIOMARKERS)) # placeholder
# # pheno_to_label_dict = dict(zip(BIOMARKERS, [p.replace('_',' ') for p in BIOMARKERS]))
# # pheno_to_label_dict['c_reactive_protein'] = 'C-reactive protein'
# -

DNAX_DATA_DIR = '/mnt/project/saige_pipeline/data'
LOCAL_WD = '/well/lindgren-ukbb/projects/ukbb-11867/nbaya/ukb_wes_450k_gwas'
LOCAL_DATA_DIR = f'{LOCAL_WD}/data'
LOCAL_PLOTS_DIR = f'{LOCAL_WD}/plots'

GWAS_DATA_DIR = LOCAL_DATA_DIR
GWAS_PLOTS_DIR = LOCAL_PLOTS_DIR 

# BIOMARKERS = pd.read_csv(
#    f'{DATA_DIR}/phenotypes/phenotype_list.qced_biomarkers.txt',
#    sep='\t',
#    names=['phenotype']
# ).phenotype.tolist()


def get_phenotype_list(phenotype_group):
    df = pd.read_csv(
        f'{GWAS_DATA_DIR}/phenotypes/phenotype_list.{phenotype_group}.txt',
        header=None
    )
    return df.values.flatten()

  
def get_obesity_phenotype_list():
    """Used for compatibility of old scripts
    """
    return get_phenotype_list(phenotype_group='obesity')


def get_saige_results_path(pheno, chrom, phenotype_group='qced_biomarkers',
                           pop='eur', sex='both_sexes', assoc='variant_test', dataset='wes'):
    dir_path = f'{GWAS_DATA_DIR}/02_saige_{assoc}{"-imputed_v3" if dataset=="imputed" else ""}/{phenotype_group}/{pop}/{sex}'
    fname = f'saige_{assoc}.{pheno}-{pop}-{sex}{"-imputed_v3" if dataset=="imputed" else ""}.chr{chrom}.tsv.gz'
    return f'{dir_path}/{fname}'


def get_old_consequence_path(chrom):
    return f'{GWAS_DATA_DIR}/annotations/old-wes_450k/ukb_wes_450k.qced.chr{chrom}.worst_csq_by_gene_canonical.txt.gz'


def get_consequence_annotation_path(chrom, version):
    version_major = int(str(version).split(".")[0])
    if chrom is None:
        # Combined across all chroms, variants grouped by variant ID and consequence category then deduplicated (combining gene IDs as needed)
        path_dir = f'{GWAS_DATA_DIR}/annotations/wes_450k'
        chrom_str = ""
    else:
        path_dir = f'/well/lindgren/barney/brava_annotation/data/vep/annotated/v{version_major}'
        chrom_str = f'.chr{chrom}'
        
    if version_major>=5:
        qc_version = 'july.qced'
    else:
        qc_version = 'qced'
        
    if version_major==6:
        file_type='tsv'
    else:
        file_type='txt'
        
    return f'{path_dir}/ukb_wes_450k.{qc_version}.brava.v{version}{chrom_str}.worst_csq_by_gene_canonical.{file_type}.gz'


def get_consequence_file(chrom, version):
    csq = pd.read_csv(
        get_consequence_annotation_path(chrom=chrom, version=version), 
        compression='gzip', 
        delim_whitespace=True,
    )
    
    csq = csq.rename(columns={'varid': 'MarkerID'})
    
    if str(version)=='5':
        raise ValueError('WARNING: Do not use v5 or v5.1 annotations (use v5.2+ instead)!!!')
#         marker_id_split = csq['MarkerID'].str.split(':', expand=True)
#         # Swap ref/alt alleles in the variant ID
#         csq['MarkerID'] = marker_id_split[0]+':'+marker_id_split[1]+':'+marker_id_split[3]+':'+marker_id_split[2]
    
    return csq

  
def get_all_consequence_files(version):
    csq_df_list = [get_consequence_file(chrom=chrom, version=version)
                   for chrom in list(range(1, 23))+['X']]
    csq = pd.concat(csq_df_list, axis=0)
    return csq


def get_unique_variant_genes_and_consequences_path(version):
    return get_consequence_annotation_path(chrom=None, version=version)


def write_unique_variant_genes_and_consequences(version, all_csq_files=None, overwrite=False):
    def concat_genes(x):
        return '/'.join(x)
    
    print(f'Writing unique-variant genes and consequences (using v{version} annotations)')

    if all_csq_files is None:
        all_csq_files = get_all_consequence_files(version=version)
    all_csq_files = all_csq_files.drop(columns=['locus', 'alleles'])
    
    # Fill missing gene symbols to avoid errors when joining missing
    all_csq_files.loc[all_csq_files['gene_symbol'].isna(), 'gene_symbol'] = 'missing'
    
    # Group variants by variant ID and consequence category then deduplicate (combining gene IDs if variant is duplicated)
    csq_grouped = all_csq_files.groupby(
        ['MarkerID', 'consequence_category']).agg((concat_genes))
    
    # Reset index to set the fields used for grouping to be column fields instead of indices
    csq_grouped = csq_grouped.reset_index()

    path = get_unique_variant_genes_and_consequences_path(version=version)
    if not isfile(path) or overwrite:
        csq_grouped.to_csv(
            path,
            compression='gzip',
            sep='\t',
            index=False
        )
        print(f'Unique-variant genes and consequences (v{version}) written to:\n{path}')

        return csq_grouped
    else:
        print(f'File already exists! (Set overwrite=True to overwrite)\n{path}')
        

def get_unique_variant_genes_and_consequences(version=CURRENT_CONSEQUENCE_VERSION):
    if str(version) in {'5','5.1','5.2'}:
        raise ValueError('WARNING: Do not use v5.x annotations (use v6+ non-unique version instead)!!!')
    print(f'Getting v{version} consequence annotations')
    csq_path = get_unique_variant_genes_and_consequences_path(version=version)
    
    if not isfile(csq_path):
        write_unique_variant_genes_and_consequences(version=version)
        
    csq_grouped = pd.read_csv(
        csq_path,
        sep='\t',
        compression='gzip'
    )

    return csq_grouped


def get_variant_genes_and_consequences(version=CURRENT_CONSEQUENCE_VERSION):
    '''
    NOTE: This differs from previous versions in that there are "duplicate" rows for 
    some variants, because they may lie in multiple genes/transcripts.
    '''
    path = get_consequence_annotation_path(chrom=None, version=version)
    last_modified_time = datetime.fromtimestamp(getmtime(path)).strftime('%Y-%m-%d')
    print(f'Getting v{version} consequence annotations (Last modified: {last_modified_time})')
    print(f'*  NOTE: Multiple rows for variants which overlap multiple genes')
    csq = pd.read_csv(
        path,
        sep='\t',
        compression='gzip'
    )
    
    if version=='6':
        csq = csq.rename(columns={
            'rsid':'MarkerID',
            'worst_csq_by_gene_canonical.gene_id': 'gene_id',
            'worst_csq_by_gene_canonical.gene_symbol': 'gene_symbol',
            'worst_csq_by_gene_canonical.most_severe_consequence': 'most_severe_csq', 
            'annotation': 'consequence_category'
        })
    
    return csq


def get_gene_id_to_symbol_map(df_dict, version=6):
    gene_id_to_symbol_map = df_dict['csq'][['gene_id','gene_symbol']].drop_duplicates()
    
    # Remove rows with missing values
    gene_id_to_symbol_map = gene_id_to_symbol_map[gene_id_to_symbol_map.count(axis=1)==2]
    
    if version < 6:
        # Get mapping of gene id to gene symbol
        # NOTE: When using df_dict['csq'] we need to account for variants which are in multiple genes
        def join_gene_id_and_symbol(row):
            joined = [f'{i};{s}' for i,s in zip(row['gene_id'].split('/'), row['gene_symbol'].split('/'))]
            return '/'.join(joined)

        gene_id_to_symbol_map['joined_id_and_symbol'] = gene_id_to_symbol_map.apply(join_gene_id_and_symbol, axis=1)

        joined_gene_id_and_symbol = gene_id_to_symbol_map.assign(joined_id_and_symbol_explode=gene_id_to_symbol_map['joined_id_and_symbol'].str.split('/')).explode('joined_id_and_symbol_explode')
        joined_gene_id_and_symbol['gene_id'] = joined_gene_id_and_symbol['joined_id_and_symbol_explode'].str.split(';', expand=True)[0]
        joined_gene_id_and_symbol['gene_symbol'] = joined_gene_id_and_symbol['joined_id_and_symbol_explode'].str.split(';', expand=True)[1]

        gene_id_to_symbol_map = joined_gene_id_and_symbol[['gene_id','gene_symbol']].drop_duplicates()
        
    gene_id_to_symbol_dict = dict(zip(*gene_id_to_symbol_map.T.values))
    
    return gene_id_to_symbol_dict

  
def get_regenie_path(pheno, chrom, test='additive', outlier_type=None):
    pheno = pheno.replace('_', '-')
    if outlier_type is None:
        return f'{GWAS_DATA_DIR}/regenie/qced_biomarkers/eur/both_sexes/step2_{test}_chr{chrom}_{pheno}_qced.regenie.gz'


def standard_variant_processing(df_dict, df, dataset, remove_ur=False):
    df = df.rename(columns={
        'POS':'position',
        'CHR':'chrom',
        'BETA':'beta',
        'SE':'se',
        'p.value':'pval'
    })
    df['chisq'] = (df['beta']/df['se'])**2
    df['pval'] = df['pval'].astype(float)
    df['nlog10pval'] = -1*np.log10(df['pval'])
    df.loc[df['pval']==0, 'nlog10pval'] = 315 # Use as placeholder for "very significant"
    
    df['maf'] = pd.concat([df['AF_Allele2'], 1-df['AF_Allele2']], axis=1).min(axis=1)
    if 'N_case' in df.columns:
        df['mac'] = pd.concat([df['AC_Allele2'], 2*(df['N_case']+df['N_ctrl'])-df['AC_Allele2']], axis=1).min(axis=1)
    
    # Only relevant for "all test"
    if remove_ur:
        is_ur_chrom = df["chrom"]=="UR"
        print(f'Removing lines with chrom = "UR" (only relevant for "saige_all_test"): {(is_ur_chrom).sum()}')
        df = df[~is_ur_chrom]
    
    df.loc[df['position']!='UR', 'position'] = df.loc[df['position']!='UR', 'position'].astype(int)
    
    if dataset=='wes':
        csq = df_dict['csq'][['MarkerID','gene_id','gene_symbol','consequence_category']]
        df = df.merge(csq, on='MarkerID', how='left')
    elif dataset=='imputed':
        print(f'Warning: imputed data may contain duplicate MarkerIDs, skipping merging with ANNOVAR')
#         csq = df_dict['annovar'][['MarkerID','region_type','gene_symbol','dist_to_gene','csq']]
    
    return df

  
def get_sample_size(n_sample_col):
    n_samples = n_sample_col.unique()
    assert len(n_samples)==1
    return n_samples[0]

  
def read_saige_gwas(df_dict, pheno, phenotype_group, pop='eur', sex='both_sexes', assoc='variant_test', dataset='wes'):

    df_list = []

    for chrom in list(range(1, 23))+['X']:
        path = get_saige_results_path(
            phenotype_group=phenotype_group,
            pheno=pheno,
            sex=sex,
            pop=pop,
            assoc=assoc,
            chrom=chrom,
            dataset=dataset
        )
        if not isfile(path):
            print(f'No file for {pheno} pop={pop} sex={sex} chr{chrom}:\n{path}')
            continue
        df_tmp = pd.read_csv(path, sep='\t', compression='gzip')
        if (chrom == 'X') and (dataset == 'wes'):
            # Needed for merging with consequence annotation file
            df_tmp['MarkerID'] = 'chr'+df_tmp[['CHR', 'POS', 'Allele2', 'Allele1']
                                              ].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)
        df_list.append(df_tmp)
    df = pd.concat(df_list, axis=0)
    df = standard_variant_processing(
        df_dict=df_dict,
        df=df,
        dataset=dataset
    )

    return df

  
def read_saige_gene_assoc(df_dict, pheno, phenotype_group, pop='eur', sex='both_sexes', assoc='set_test'):

    df_list = []

    for chrom in list(range(1, 23))+['X']:
        path = get_saige_results_path(
            pheno=pheno,
            phenotype_group=phenotype_group,
            chrom=chrom,
            sex=sex,
            assoc=assoc
        )
        if not isfile(path):
            print(f'No file for {pheno} sex={sex} chr{chrom}: {path}')
            continue
        df_tmp = pd.read_csv(path, sep='\t')
        df_tmp['chrom'] = str(chrom)
        for pval_field in ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']:
            df_tmp[f'nlog10{pval_field}'] = -np.log10(df_tmp[pval_field])
        df_list.append(df_tmp)
    df = pd.concat(df_list, axis=0)

    gene_id_to_symbol_df = df_dict['csq'][[
        'gene_symbol', 'gene_id']].drop_duplicates(subset='gene_id')
    df = df.merge(gene_id_to_symbol_df, left_on='Region',
                  right_on='gene_id', how='left')

    return df


def read_regenie_results(pheno):
    def get_path(chrom): return get_regenie_path(pheno=pheno, chrom=chrom)
    df_list = [pd.read_csv(get_path(chrom), compression='gzip', sep=' ')
               for chrom in range(1, 23)]
    df = pd.concat(df_list, axis=0)

    df = df.rename(columns={
        'CHROM': 'chrom',
        'GENPOS': 'position',
        'ID': 'MarkerID',
        'LOG10P': 'nlog10pval',
        'ALLELE0': 'Allele1',
        'ALLELE1': 'Allele2',
        'A1FREQ': 'AF_Allele2',
    })

    return df


def read_genebass_results(pheno, assoc='variant'):
    genebass_suffix = genebass_file_dict[pheno]
    if assoc == 'gene':
        gb = pd.read_csv(
            f'{GWAS_DATA_DIR}/genebass/gene-burden-results-exomes_pLoF_continuous-{genebass_suffix}.csv')
    elif assoc == 'variant':
        gb = pd.read_csv(
            f'{GWAS_DATA_DIR}/genebass/single-variant-associations-exomes_continuous-{genebass_suffix}.csv')
        gb['CHR'] = gb['Variant ID'].str.split('-', expand=True)[0]
        gb['position'] = gb['Variant ID'].str.split('-', expand=True)[1]
        gb['Allele1'] = gb['Variant ID'].str.split('-', expand=True)[2]
        gb['Allele2'] = gb['Variant ID'].str.split('-', expand=True)[3]

    return gb

  
def is_coding_nonsynonymous(df):
    return (
        df.consequence_category.notna()
        & ~df.consequence_category.isin({'non_coding', 'synonymous'})
    )


def get_imputed_v3_hwe_pval_path(chrom):
    return f'{GWAS_DATA_DIR}/annotations/imputed_v3/ukb_impv3.subset_to_wes_450k_qc_pass_eur.chr{chrom}.hardy.gz'


def get_imputed_v3_hwe_pval(chrom_list=list(range(1, 23))+['X']):
    df_list = []
    for chrom in chrom_list:
        df_chrom = pd.read_csv(
            get_imputed_v3_hwe_pval_path(chrom),
            delim_whitespace=True,
            compression='gzip'
        )
        df_list.append(df_chrom)

    df = pd.concat(df_list, axis=0)
    df = df.rename(columns={
        '#CHROM': 'chrom',
        'P': 'hwe_pval',
        'ID': 'MarkerID'
    })

    return df


def get_heterogeneity_se(beta1, beta2, se1, se2):
    if any(np.isnan(beta1)) or any(np.isnan(beta2)):
        raise ValueError('Warning: There exist undefined values in the beta values')
    spearman_r, pval = stats.spearmanr(beta1, beta2)
    se_diff = np.sqrt(np.square(se1) + np.square(se2) - 2*spearman_r*se1*se2)
    return se_diff


def get_heterogeneity_tstat(beta1, beta2, se1, se2):
    beta_diff = beta1-beta2
    se_diff = get_heterogeneity_se(beta1, beta2, se1, se2)

    return beta_diff/se_diff


def get_heterogeneity_gwas(gwas1, gwas2, suffixes=['_f', '_m'], heterogeneity_suffix='_sexdiff', 
    beta_field='beta', se_field='se', pval_field='pval',
    merge_on=['chrom', 'position', 'MarkerID', 'Allele1', 'Allele2',
             'gene_id', 'gene_symbol', 'consequence_category']):
    merge = gwas1.merge(
        gwas2,
        on=merge_on,
        suffixes=suffixes
    )
    
    merge[f'{beta_field}{heterogeneity_suffix}'] = merge[f'{beta_field}{suffixes[0]}']-merge[f'{beta_field}{suffixes[1]}']

    merge[f'tstat{heterogeneity_suffix}'] = get_heterogeneity_tstat(
        beta1=merge[f'{beta_field}{suffixes[0]}'],
        beta2=merge[f'{beta_field}{suffixes[1]}'],
        se1=merge[f'{se_field}{suffixes[0]}'],
        se2=merge[f'{se_field}{suffixes[1]}']
    )

    merge[f'{pval_field}{heterogeneity_suffix}'] = 2 * \
        stats.norm.cdf(-merge[f'tstat{heterogeneity_suffix}'].abs())
    merge[f'nlog10pval{heterogeneity_suffix}'] = - \
        np.log10(merge[f'{pval_field}{heterogeneity_suffix}'])

    return merge

  
def get_brava_group_test_path(chrom, phenotype_group='obesity', pop='EUR', pheno='bmi'):
    if pheno=='body_mass_index_bmi':
        pheno='BMI'
    return f'{GWAS_DATA_DIR}/02_saige_all_test-brava/{phenotype_group}/{pop.lower()}/both_sexes/chr{chrom}_{pheno}_{pop.upper()}.txt.gz'

  
def get_brava_variant_test_path(chrom, phenotype_group='obesity', pop='EUR', pheno='bmi'):
    if pheno=='body_mass_index_bmi':
        pheno='BMI'
    return f'{GWAS_DATA_DIR}/02_saige_all_test-brava/{phenotype_group}/{pop.lower()}/both_sexes/chr{chrom}_{pheno}_{pop.upper()}.txt.singleAssoc.txt.gz'

  
def get_all_test_variant_assoc_path(pheno, chrom, phenotype_group, pop, sex):
    return f'{GWAS_DATA_DIR}/02_saige_all_test/{phenotype_group}/{pop}/{sex}/saige_all_test.{pheno}-{pop}-{sex}.chr{chrom}.tsv.singleAssoc.txt.gz'
    
    
def load_saige_all_test_variant_assoc(df_dict, gwas_id, phenotype_group='obesity'):
    all_test_variant_assoc_id = f'{gwas_id}-all_test_variants'
    if all_test_variant_assoc_id not in df_dict:
        pheno, pop, sex = gwas_id.split('-')

        group_test_df_list = []
        variant_test_df_list = []
        for chrom in list(range(1,23))+['X']:
            path = get_all_test_variant_assoc_path(
                pheno=pheno,
                phenotype_group=phenotype_group,
                pop=pop,
                sex=sex,
                chrom=chrom
            )
            try:
                variant_test_df_list.append(pd.read_csv(
                    path,
                    sep='\t',
                    compression='gzip'
                ))
            except: 
                print(f'File does not exist: {path}')


        variant_assoc = pd.concat(variant_test_df_list, axis=0)
        variant_assoc = standard_variant_processing(
            df_dict=df_dict, 
            df=variant_assoc,
            dataset='wes'
        )
        
        df_dict[all_test_variant_assoc_id] = variant_assoc
    
    return df_dict[all_test_variant_assoc_id]

  
def load_finemapping_results(df_dict, gwas_id):
    master_id = f'{gwas_id}-master'
    imputed_gwas_id = f'{gwas_id}-imputed'
    if master_id not in df_dict:
        ldstore_master = pd.read_csv(
            get_ldstore_master_path(imputed_gwas_id),
            sep=';',
        )

        ldstore_master['chrom'] = ldstore_master['bgen'].str.split('chr', expand=True)[1].str.replace('.bgen','')
        ldstore_master['locus_id'] = [z.split('/')[-1].strip('.z') for z in ldstore_master['z'].tolist()]

        df_dict[master_id] = ldstore_master

    ldstore_master = df_dict[master_id]

    print(f'Loci in {master_id}: {ldstore_master.shape[0]}')
    
    return df_dict, df_dict[master_id]


def load_imputed_results(df_dict, phenotype_group, gwas_id, force_read=False):
    from association_processing.finemap_utils import read_finemap_output # put inside this function to avoid circular import
    
    pheno, pop, sex = gwas_id.split('-')
    imputed_gwas_id = f'{gwas_id}-imputed'
    print(imputed_gwas_id)
    if (imputed_gwas_id not in df_dict) or force_read:
        df_dict[imputed_gwas_id] = read_saige_gwas(
            df_dict=df_dict, 
            pheno=pheno, 
            phenotype_group=phenotype_group, 
            pop=pop, 
            sex=sex,
            dataset='imputed'
        )

    if 'prob' not in df_dict[imputed_gwas_id].columns:
        finemapped_loci = read_finemap_output(df_dict=df_dict, locus_id=imputed_gwas_id)
        finemapped_loci = finemapped_loci.rename(columns={'rsid':'MarkerID'})
        finemapped_loci = finemapped_loci[['MarkerID','prob','lead_variant']]
        df_dict[imputed_gwas_id] = df_dict[imputed_gwas_id].copy().merge(finemapped_loci, on='MarkerID', how='left')
    
#     if gene_symbol_in_cols:
#         df_dict[imputed_gwas_id] = df_dict[imputed_gwas_id].drop(columns=['gene_symbol'])
    
#     if 'hgnc_symbol' not in df_dict[imputed_gwas_id].columns:
#         df_dict[imputed_gwas_id] = df_dict[imputed_gwas_id].merge(df_dict['imputed_v2g'], on='MarkerID', how='left')
        
    return df_dict, df_dict[imputed_gwas_id]
    
    
def load_wes_results(df_dict, gwas_id, phenotype_group):
    pheno, pop, sex = gwas_id.split('-')
    wes_gwas_id = f'{gwas_id}-wes'
    print(wes_gwas_id)
    if wes_gwas_id not in df_dict:
        df_dict[wes_gwas_id] = read_saige_gwas(
                df_dict=df_dict, 
                pheno=pheno, 
                phenotype_group=phenotype_group, 
                pop=pop, 
                sex=sex,
                assoc='variant_test'
            )

    return df_dict, df_dict[wes_gwas_id]

  
def prepare_sexdiff(df_dict, gwas_id, phenotype_group, group_test='pLoF;damaging_missense', max_maf=0.01):
    pheno, pop, sex = gwas_id.split('-')
    from association_processing.finemap_utils import read_finemap_output # put inside this function to avoid circular import
     
    imputed_gwas_id = f'{gwas_id}-imputed'
    wes_gwas_id = f'{gwas_id}-wes'
    set_gwas_id = f'{gwas_id}-set'
    
    if imputed_gwas_id not in df_dict:
        for tmp_sex in ['male','female']:
            tmp_imputed_gwas_id = imputed_gwas_id.replace('-sexdiff-',f'-{tmp_sex}-')
            print(tmp_imputed_gwas_id)
            if tmp_imputed_gwas_id not in df_dict:
                df_dict[tmp_imputed_gwas_id] = read_saige_gwas(
                    df_dict=df_dict, 
                    pheno=pheno, 
                    phenotype_group=phenotype_group, 
                    pop=pop, 
                    sex=tmp_sex,
                    dataset='imputed'
                )
        male = df_dict[imputed_gwas_id.replace('-sexdiff-',f'-male-')]
        female = df_dict[imputed_gwas_id.replace('-sexdiff-',f'-female-')]

        df_dict[imputed_gwas_id] = get_heterogeneity_gwas(male, female)

    finemapped_loci = read_finemap_output(df_dict=df_dict, locus_id=imputed_gwas_id)
    finemapped_loci = finemapped_loci.rename(columns={'rsid':'MarkerID'})
    finemapped_loci = finemapped_loci[['MarkerID','prob']]

    if 'prob' not in df_dict[imputed_gwas_id].columns:
        df_dict[imputed_gwas_id] = df_dict[imputed_gwas_id].copy().merge(finemapped_loci, on='MarkerID', how='left')
        
    if 'hgnc_symbol' not in df_dict[imputed_gwas_id].columns:
        df_dict[imputed_gwas_id] = df_dict[imputed_gwas_id].merge(df_dict['imputed_v2g'], on='MarkerID', how='left')

    imputed_gwas = df_dict[imputed_gwas_id]

    if wes_gwas_id not in df_dict:
        for tmp_sex in ['male','female']:
            tmp_wes_gwas_id = wes_gwas_id.replace('-sexdiff-',f'-{tmp_sex}-')
            print(tmp_wes_gwas_id)
            if tmp_wes_gwas_id not in df_dict:
                df_dict[tmp_wes_gwas_id] = read_saige_gwas(
                    df_dict=df_dict, 
                    pheno=pheno, 
                    phenotype_group=phenotype_group, 
                    pop=pop, 
                    sex=tmp_sex,
                )
        male = df_dict[wes_gwas_id.replace('-sexdiff-',f'-male-')]
        female = df_dict[wes_gwas_id.replace('-sexdiff-',f'-female-')]

        df_dict[wes_gwas_id] = get_heterogeneity_gwas(male, female)

    # Update fields for plotting
    for tmp_gwas_id in [imputed_gwas_id, wes_gwas_id]: 
        df = df_dict[tmp_gwas_id]
        for field in ['AC_Allele2','AF_Allele2']:
            df[field] = df[[f'{field}_m',f'{field}_f']].mean(axis=1)
        df['maf'] = pd.concat([df['AF_Allele2'], 1-df['AF_Allele2']], axis=1).min(axis=1)
        df['beta'] = (df['beta_f']-df['beta_m']).abs()
        df['pval'] = df['pval_sexdiff']
        df['nlog10pval'] = df['nlog10pval_sexdiff']
        df_dict[tmp_gwas_id] = df

    if set_gwas_id not in df_dict:
        for tmp_sex in ['male','female']:
            tmp_set_gwas_id = set_gwas_id.replace('-sexdiff-',f'-{tmp_sex}-')
            print(tmp_set_gwas_id)
            if tmp_set_gwas_id not in df_dict:
                df_dict[tmp_set_gwas_id] = read_saige_gene_assoc(
                    df_dict=df_dict, 
                    pheno=pheno, 
                    phenotype_group=phenotype_group, 
                    pop=pop, 
                    sex=tmp_sex,
                )
        male = df_dict[set_gwas_id.replace('-sexdiff-',f'-male-')]
        female = df_dict[set_gwas_id.replace('-sexdiff-',f'-female-')]
        female = female[female['Group']==group_test]
        male = male[male['Group']==group_test]

        df_dict[set_gwas_id] = get_heterogeneity_gwas(
            gwas1=female, 
            gwas2=male,  
            heterogeneity_suffix='_sexdiff', 
            beta_field='BETA_Burden', 
            se_field='SE_Burden',
            merge_on=['Region', 'Group', 'max_MAF', 'gene_symbol', 'gene_id','CHR'])

    return df_dict

  
def get_significant_variant_test_results_path(phenotype_group, suffix):
    return f'{GWAS_DATA_DIR}/significant_results/ukb_wes_450k.{phenotype_group}.varianttest.{suffix}.tsv.gz'

  
def write_significant_variant_test_results(df, phenotype_group, suffix, overwrite=False):
    path = get_significant_variant_test_results_path(phenotype_group=phenotype_group, suffix=suffix)
    
    if not isfile(path) or overwrite:
        print(f'Writing variant test results:\n{path}')
        print(f'Number of rows: {df.shape[0]}')
        df.to_csv(
            path,
            sep='\t',
            compression='gzip',
            index=False
        )
    else:
        print(f'File exists! (use overwrite=True to overwrite):\n{path}')
    
    
def get_significant_set_test_results_path(phenotype_group, group_test, max_maf, suffix=None):
    if suffix is None:
        if group_test != 'Cauchy':
            max_maf_str = '.maxmaf_{max_maf}'
        else:
            max_maf_str = ''
        suffix = f'{group_test.replace(";","-")}{max_maf_str}'

    return f'{GWAS_DATA_DIR}/significant_results/ukb_wes_450k.{phenotype_group}.grouptest.{suffix}.tsv.gz'

  
def write_significant_set_test_results(df, phenotype_group, group_test, max_maf, suffix=None, overwrite=False):
    path = get_significant_set_test_results_path(
        phenotype_group=phenotype_group, 
        group_test=group_test, 
        max_maf=max_maf,
        suffix=suffix
    )
    
    if suffix is None:
        if (group_test=='Cauchy') and ('MAF_aggregated' in df.columns):
            df = df.drop(columns='MAF_aggregated')
        
    if ('Region' in df.columns) & ('gene_id' in df.columns):
        df = df.drop(columns='gene_id')
        
    if not isfile(path) or overwrite:
        print(f'Writing gene test results:\n{path}')
        print(f'Number of rows: {df.shape[0]}')
        df.to_csv(
            path,
            sep='\t',
            compression='gzip',
            index=False
        )
    else:
        print(f'File exists! (use overwrite=True to overwrite):\n{path}')
 

def remove_noncoding_synonymous(df):
    '''NOTE: Only use when working with data with incomplete annotations. 
    This is not a perfect solution.
    '''
    print('Removing non-coding or synonymous variants')
    n_initial = df.shape[0]
    print(f'* Initial variant count: {n_initial}')
    n_missing = df["consequence_category"].isna().sum()
    print(f'* Variants with missing consequence category: {n_missing} ({100*n_missing/n_initial:.2f}%)')
    is_noncoding_or_synonymous = df['consequence_category'].isin(['synonymous', 'non_coding'])
    n_removed = is_noncoding_or_synonymous.sum()
    print(f'* Variants removed (non-coding or synonymous): {n_removed} ({100*n_removed/n_initial:.2f}%)')
    df = df[~is_noncoding_or_synonymous]
    

    print(f'* Post-filter variant count: {df.shape[0]}')
    
    return df

  
def annotate_with_significance_field(df, genomewide_sig, nominally_sig, remove_nonsig=False, pval_field='pval'):
    print("Annotating with 'significance' field")
    
    # Add 'significance' field

    print(f"* Variants missing '{pval_field}' values: {df[pval_field].isna().sum()}")
    is_genomewide_sig = df[pval_field] < genomewide_sig
    # NOTE: We exclude genome-wide significant variants from nominally significant variants
    is_nominally_sig = (
        (df[pval_field] < nominally_sig)
        & (~is_genomewide_sig)
    )
    df.loc[is_genomewide_sig, 'significance'] = 'genomewide_significant'
    df.loc[is_nominally_sig, 'significance'] = 'nominally_significant'
    
    print(f'* Significant variants:')
    print(f'  - Genome-wide significant (P<{genomewide_sig}): {is_genomewide_sig.sum()}')
    print(f'  - Nominally significant (P>={genomewide_sig} & P<{nominally_sig}): {is_nominally_sig.sum()}')
    
    if remove_nonsig:
        print(f'Removing non-significant variants (remove_nonsig=True)')
        print(f'* Initial variant count: {df.shape[0]}')
        df = df[df['significance'].notna()]
        print(f'* Post-filter variant count: {df.shape[0]}')
    
    return df


def print_summary(df, phenos, assoc_test='variant'):
    assert assoc_test in {'variant','gene'}
    if assoc_test=='variant':
        drop_duplicates_subset='MarkerID'
    elif assoc_test=='gene':
        drop_duplicates_subset='Region'
    n_total = df.shape[0]
    n_unique = df.drop_duplicates(subset=drop_duplicates_subset).shape[0]
    print(f'Total rows: {n_total}')
    print(f'* Unique {assoc_test}s: {n_unique}')
    
    # Per-sex breakdown
    print('\nPer sex:')
    for sex in ['both_sexes','female','male']:
        n_per_sex = (df["sex"]==sex).sum()
        print(f'* {sex}: {n_per_sex} ({100*n_per_sex/n_total:.2f}%)')

    # Per-phenotype breakdown
    print('\nPer phenotype:')
    for pheno in phenos:
        n_per_pheno = (df["pheno"]==pheno).sum()
        print(f'* {pheno}: {n_per_pheno} ({100*n_per_pheno/n_total:.2f}%)')

    # Per-chrom breakdown
    print('\nPer chromosome:')
    for chrom in list(range(1,23))+['X']:
        n_per_chrom = (df["chrom"].astype(str)==str(chrom)).sum()
        print(f'* {chrom}: {n_per_chrom} ({100*n_per_chrom/n_total:.2f}%)')

    if 'consequence_category' in df.columns:
        # Per-consequence breakdown
        print('\nPer consequence:')
        for csq in df['consequence_category'].unique():
            if str(csq)=='nan':
                n_per_csq = (df["consequence_category"].isna()).sum()
            else:
                n_per_csq = (df["consequence_category"]==csq).sum()
            print(f'* {csq}: {n_per_csq} ({100*n_per_csq/n_total:.2f}%)')
            
    if 'Group' in df.columns:
        # Per-group breakdown
        print('\nPer group:')
        for g in df['Group'].unique():
            print(f'* {g}: {(df["Group"]==g).sum()}')
            
    if 'max_MAF' in df.columns:
        # Per-max-MAF breakdown
        print('\nPer max MAF:')
        for max_maf in [0.01,0.001,0.0001]:
            print(f'* {max_maf}: {(df["max_MAF"]==max_maf).sum()}')
            
    if 'maf' in df.columns:
        print('\nPer MAF bin:')
        maf_thresholds = [0, 1e-6, 1e-4, 1e-2, 0.5]
        maf_bins = zip(maf_thresholds[:-1], maf_thresholds[1:])
        for maf_lower, maf_upper in maf_bins:
            n_bin = ((df["maf"]>maf_lower)&(df["maf"]<=maf_upper)).sum()
            print(f'* ({maf_lower}, {maf_upper}]: {n_bin} ({100*n_bin/n_total:.2f}%)')
      
    pip_col_names = [c for c in ['PIP','prob'] if c in df.columns]
    if len(pip_col_names)>0:
        print('\nPer PIP bin:')
        pip_col_name = pip_col_names[0]
        pip_thresholds = [0, 0.1, 0.5, 0.9, 0.99, 1]
        pip_bins = zip(pip_thresholds[:-1], pip_thresholds[1:])
        for pip_lower, pip_upper in pip_bins:
            n_bin = ((df[pip_col_name]>pip_lower)&(df[pip_col_name]<=pip_upper)).sum()
            print(f'* ({pip_lower}, {pip_upper}]: {n_bin} ({100*n_bin/n_total:.2f}%)')
            

# def get_duncan_transformed_burden_test_betas_maf(df_dict, gwas_id, phenotype_group, csq_categories_list, max_maf_list = [0.01,0.001,0.0001]):
    
#     all_test_rare_variants = load_saige_all_test_variant_assoc(
#         df_dict=df_dict,
#         gwas_id=gwas_id,
#         phenotype_group=phenotype_group
#     )
#     all_test_rare_variants['gene_id_split'] = all_test_rare_variants['gene_id'].str.split('/')
#     all_test_rare_variants = all_test_rare_variants.explode('gene_id_split')

#     # filter to variants:
#     # - MAF <= max_maf
#     # - consequence = "pLoF" or "damaging_missense"

#     grouped_df_list = []
    
#     for max_maf in max_maf_list:
#         for csq_categories in csq_categories_list:

#             group = ';'.join(csq_categories)

#             subset = all_test_rare_variants[
#                 (all_test_rare_variants.maf<=max_maf)
#                 & (all_test_rare_variants.consequence_category.isin(csq_categories))
#             ]

#             # Sum across all variants in gene
#             subset_grouped = subset.groupby(['gene_id_split']).sum()

#             subset_grouped['BETA_burden_duncan'] = subset_grouped['Tstat']/subset_grouped['var']
#             subset_grouped['Group'] = group
#             subset_grouped['max_MAF'] = max_maf
#             subset_grouped = subset_grouped.reset_index()

#             grouped_df_list.append(subset_grouped)
        
#     grouped = pd.concat(grouped_df_list, axis=0)
    
#     grouped = grouped.rename(
#         columns={
#             'gene_id': 'old_gene_id',
#             'gene_id_split':'gene_id',
#             'maf': 'burden_maf'
#         }
#     )
    
#     return grouped


def get_duncan_transformed_burden_test_betas_maf(df_dict, gwas_id, phenotype_group, csq_categories_list, max_maf_list = [0.01,0.001,0.0001]):
    
    all_test_rare_variants = load_saige_all_test_variant_assoc(
        df_dict=df_dict,
        gwas_id=gwas_id,
        phenotype_group=phenotype_group
    )
    
    if any(all_test_rare_variants['gene_id'].str.contains('/')):
        all_test_rare_variants['gene_id_split'] = all_test_rare_variants['gene_id'].str.split('/')
        all_test_rare_variants = all_test_rare_variants.explode('gene_id_split')
        all_test_rare_variants = all_test_rare_variants.rename(columns={
            'gene_id': 'old_gene_id',
            'gene_id_split':'gene_id',
        })

    # Extract and process aggregated ultra-rare ("UR") variant association results
    ur = all_test_rare_variants[all_test_rare_variants['chrom']=='UR']
    ur['gene_id'] = ur['MarkerID'].str.split(':', expand=True)[0]
    ur['Group'] = ur['MarkerID'].str.split(':', expand=True)[1]
    ur['max_MAF'] = ur['MarkerID'].str.split(':', expand=True)[2].astype(float)

    # For each combination of max_MAF and consequence group, group variants and sum 

    grouped_df_list = []
    
    for max_maf in max_maf_list:
        for csq_categories in csq_categories_list:
            
            subset = all_test_rare_variants[
                (all_test_rare_variants.maf<=max_maf)
                & (all_test_rare_variants.consequence_category.isin(csq_categories))
            ][['gene_id','maf','Tstat','var']] # for efficiency, only extract necessary columns
            
            group = ';'.join(csq_categories)
            
            ur_group = ur[
                (ur['Group']==group)
                & (ur['max_MAF']==max_maf)
            ][['gene_id','maf','Tstat','var']]
            
            # Append ur_group to the end
            subset = pd.concat([subset, ur_group], axis=0)

            # Sum across all variants in gene
#             subset_grouped = subset.groupby(['gene_id']).sum()
            subset_grouped = subset.groupby(['gene_id']).sum() # for efficiency, only group necessary columns

#             subset_grouped_ct = subset.groupby(['gene_id']).count()
#             genes_with_multiple_variants = subset_grouped_ct[subset_grouped_ct['MarkerID']>1].index
#             subset_grouped = subset_grouped[subset_grouped.index.isin(genes_with_multiple_variants)]

            subset_grouped['BETA_burden_duncan'] = subset_grouped['Tstat']/subset_grouped['var']
            
            subset_grouped['Group'] = group
            subset_grouped['max_MAF'] = max_maf
            subset_grouped = subset_grouped.reset_index()

            grouped_df_list.append(subset_grouped)
        
    grouped = pd.concat(grouped_df_list, axis=0)
    
    grouped = grouped.rename(
        columns={
            'maf': 'burden_maf'
#             'AF_Allele2': 'burden_maf'
        }
    )
    
    return grouped


# def get_burden_beta_and_maf_fields(df_dict, gwas_id, burden_test, 
#     csq_categories=[['pLoF'],['pLoF','damaging_missense']]):
    
#     variant_assoc = load_saige_all_test_variant_assoc(
#         df_dict=df_dict,
#         gwas_id=gwas_id
#     )
    
#     grouped = get_duncan_transformed_burden_tests(wes=variant_assoc)
    
#     merged_w_grouped = burden_test.merge(
#         grouped[['Group','gene_id','maf_grouped','BETA_burden_duncan']], 
#         on=['Group','gene_id'], 
#         how='left'
#     )
    
#     return merged_w_grouped

# def deprecated_get_gene_aggregated_maf(df_dict, gwas_id, group_test='pLoF;damaging_missense', max_maf=0.01):
#     set_gwas_id = f'{gwas_id}-set'
#     print(set_gwas_id)
#     if set_gwas_id not in df_dict:
#         df_dict[set_gwas_id] = read_saige_gene_assoc(
#             df_dict=df_dict, 
#             pheno=pheno, 
#             phenotype_group=phenotype_group,
#             pop=gwas_id.split('-')[1],
#             sex=gwas_id.split('-')[2],
#         )

#     genes = df_dict[set_gwas_id]

#     if 'MAF_aggregated' not in genes.columns:
#         wes_gwas_id = f'{gwas_id}-wes'

#         # Choose set of results to use for aggregating MAF
#         max_set = genes[
#             (genes['Group']==group_test)
#             &(genes['max_MAF']==max_maf)
#         ]

#         wes = df_dict[wes_gwas_id]
#         if 'sexdiff' in gwas_id:
#             if ('N' not in wes.columns) and all([c for c in ['N_f','N_m'] if c in wes.columns]):
#                 wes['N'] = wes['N_f']+wes['N_m']
#             for field in ['MAC','Number_rare','Number_ultra_rare']:
#                 max_set[field] = max_set[f'{field}_m'] + max_set[f'{field}_f']
        
#         n_samples = get_sample_size(wes['N'])
#         max_set['MAF_aggregated'] = max_set['MAC']/((max_set['Number_rare']+max_set['Number_ultra_rare'])*2*n_samples)

#         # Annotate *all* group tests with the same aggregated MAF from the single group test we've chosen

#         df_dict[set_gwas_id] = genes.merge(max_set[['Region', 'MAF_aggregated']], on='Region')
    
#     return df_dict, df_dict[set_gwas_id]

def get_effect_size_of_minor_allele(df, eaf_col='AF_Allele2', maf_col='maf', beta_col='beta'):
    """
    eaf_col : Effect allele (allele used for genotype value in regression) frequency column
    """
    sign = (
        -1*(df[eaf_col]==(1-df[maf_col]))+
        (df[eaf_col]==(df[maf_col]))+
        (df[maf_col]==0.5) # necessary to make sign nonzero when MAF=0.5
    ) 
    
    assert all(sign!=0)
    return sign * df[beta_col]

                                                          
def get_gene_aggregated_maf_and_burden_beta(df_dict, phenotype_group, gwas_id, 
    csq_categories_list = [['pLoF'],['pLoF','damaging_missense']]):
    all_test_genes_gwas_id = f'{gwas_id}-all_test_genes'
    
    print(all_test_genes_gwas_id)
    if all_test_genes_gwas_id not in df_dict:
        pheno, pop, sex = gwas_id.split('-')
        df_dict[all_test_genes_gwas_id] = read_saige_gene_assoc(
            df_dict=df_dict, 
            pheno=pheno, 
            phenotype_group=phenotype_group,
            pop=pop,
            sex=sex,
            assoc='all_test'
        )

    genes = df_dict[all_test_genes_gwas_id]
    
    if 'burden_maf' not in genes.columns:
        
        grouped = get_duncan_transformed_burden_test_betas_maf(
            df_dict=df_dict, 
            gwas_id=gwas_id, 
            phenotype_group=phenotype_group, 
            csq_categories_list=csq_categories_list, 
            max_maf_list = genes['max_MAF'].unique()
        )
    
        genes = genes.merge(
            grouped[['Group','max_MAF','gene_id','burden_maf','BETA_burden_duncan']], 
            on=['Group','max_MAF','gene_id'], 
            how='left'
        )

    df_dict[all_test_genes_gwas_id] = genes
    
    return df_dict, df_dict[all_test_genes_gwas_id]

  
def get_annovar_path(chrom, suffix='variant'):
    assert suffix in {'variant','exonic_variant'}
    return f'{GWAS_DATA_DIR}/annotations/imputed_v3/annovar/imputed_v3.chr{chrom}.{suffix}_function'

  
def get_annovar_result(chrom, suffix):
    assert suffix in {'variant','exonic_variant'}
    
    if suffix=='variant':
        names=['region_type','gene_info']
    else:
        names=['line_number','csq variant_type','gene_info']
        
    annovar = pd.read_csv(
        get_annovar_path(chrom, suffix=suffix),
        sep='\t',
        names=names+['variant_id']
    )
    
    if suffix=='exonic_variant':
        csq_variant_type_split = annovar['csq variant_type'].str.split(' ', expand=True)
        annovar['csq'] = csq_variant_type_split[0]
        annovar['variant_type'] = csq_variant_type_split[1]
        annovar = annovar[['line_number','csq','variant_type','gene_info','variant_id']]
    
#     if extract_fields_from_variant_id:
    split_varid = annovar['variant_id'].str.split(' ', expand=True)
    for i, new_col_name in enumerate(['chrom','variant_start','variant_end','ref','alt','rsid']):
        if new_col_name in {'variant_start','variant_end'}:
            split_varid[i] = split_varid[i].astype(int)
        annovar[new_col_name] = split_varid[i]
    annovar = annovar.drop(columns='variant_id')
    
    return annovar

  
def get_merged_annovar_result(chrom, rename_cols_to_match_saige=True):
    func = get_annovar_result(chrom=chrom, suffix='variant')
    csq = get_annovar_result(chrom=chrom, suffix='exonic_variant')

    merged = func.merge(
        csq, 
        on=['chrom','variant_start','variant_end','ref','alt','rsid'], 
        how='left',
        suffixes=('','_exonic')
    )
    
    if rename_cols_to_match_saige:
        merged = merged.rename(
            columns={
                'variant_start':'position',
                'ref':'Allele1',
                'alt':'Allele2',
                'rsid':'MarkerID'
            }
        )
    
    return merged

  
def explode_gene_info(df):
#     df = df.drop_duplicates(subset=['region_type','gene_info','chrom','position','Allele1','Allele2'])

    # For all region_types with gene_info values that DO NOT include ":", split on comma and explode
    is_not_contains_colon = ~df['gene_info'].str.contains(':')
    df.loc[is_not_contains_colon, 'gene_info_split'] = df.loc[is_not_contains_colon, 'gene_info'].str.split(',')
    df = df.explode('gene_info_split')

    # If gene info contains distance, extract the distance
    is_contains_dist = df['gene_info_split'].notna() & df['gene_info_split'].str.contains('dist=')
    df.loc[is_contains_dist, 'dist_to_gene'] = df.loc[is_contains_dist, 'gene_info_split'].apply(lambda s: re.search('dist=(.*)\)', s).group(1))

    # For the remaining variants, split on ")," and explode
    is_not_split = df['gene_info_split'].isna()
    df.loc[is_not_split, 'gene_info_split'] = df.loc[is_not_split, 'gene_info'].str.split('\),')
    df = df.explode('gene_info_split')

    # Remove anything in open brackets
    df['gene_symbol'] = df['gene_info_split'].apply(lambda x: re.sub(r"\(.*", "", x))
    
    # Remove unnecessary column
    df = df.drop(columns={'gene_info_split'})
    
    return df

  
def get_ensembl_gene_and_transcript_info(build='GRCh37'):
    assert build in {'GRCh37'}
    print(f'Using build {build}')
    path = f'{GWAS_DATA_DIR}/annotations/imputed_v3/ensembl_transcripts.{build}.tsv.gz'
    ensembl = pd.read_csv(
        path,
        sep='\t',
        compression='gzip'
    )
    
    return ensembl


def get_exploded_annovar_path(version):
    version = version.strip('v')
    return f'{GWAS_DATA_DIR}/annotations/imputed_v3/imputed_v3.annovar_gene_and_csq_annotations.v{version}.tsv.gz'


def write_exploded_annovar(version):
    version = version.strip('v')
    
    if version=='1.0':
        annovar_df_list = []

        for chrom in list(range(1,23))+['X']:
            annovar_df_list.append(get_merged_annovar_result(chrom=chrom))

        annovar = pd.concat(annovar_df_list, axis=0)

        annovar_explode = explode_gene_info(annovar)
        
        annovar_explode = annovar_explode[[
            'chrom','position','variant_end','Allele1','Allele2',
            'MarkerID','region_type','gene_symbol','dist_to_gene','csq'
        ]]
        
        path = get_exploded_annovar_path(version=version)
        
        annovar_explode.to_csv(
            path,
            sep='\t',
            compression='gzip',
            index=False
        )
        print(f'Exploded ANNOVAR file written to:\n{path}')
        
        return annovar_explode
    else:
        raise ValueError(f'Version not supported: {version}')

        
def get_exploded_annovar(version='1.0'):
    print(f'Getting ANNOVAR variant-gene annotations (v{version})')
    path = get_exploded_annovar_path(version=version)
    
    if isfile(path):
        annovar_explode = pd.read_csv(
            path, 
            sep='\t',
            compression='gzip',
            dtype={
                'chrom':str,
                'dist_to_gene': float,
            },
            na_values=['NONE']
        )
    else:
        print(f'Exploded ANNOVAR file does not yet exist:\n{path}\n\nWriting new file...')
        annovar_explode = write_exploded_annovar(version=version)
    return annovar_explode


def get_hwat_results():
    # Read hWAT results
    hwat = pd.read_csv(
        f'{GWAS_DATA_DIR}/rna_seq/rna_seq_adipo_anno_count_normalised_hWAt.txt.gz',
        compression='gzip',
        delim_whitespace=True
    )

    hwat = hwat.rename(columns={'gene_id':'gene_symbol'})
    hwat = hwat.reset_index(names=['gene_id'])
    hwat = hwat[['gene_id','gene_symbol']+[c for c in hwat.columns if 'hWAT' in c]]
    # Drop gene_symbol, we will merge our results on gene_id instead
    hwat = hwat.drop(columns='gene_symbol')

    days = sorted(set([int(c.split('.')[0].strip('hWAT_D')) for c in hwat.columns if c[:6]=='hWAT_D']))

    for day in days:
        # for each day, annotate with mean and std err of mean across three replicates
        for stat, stat_fn in [('mean',pd.DataFrame.mean), ('sem',pd.DataFrame.sem)]:
            hwat[f'hWAT_D{day}_{stat}'] = stat_fn(hwat[[f'hWAT_D{day}{rep}' for rep in ['','.1','.2']]], axis=1)

        # if significantly less than 10 (95% CI for mean is less than 10)
        hwat[f'hWAT_D{day}_sig_lt_10'] = (hwat[f'hWAT_D{day}_mean']+1.96*hwat[f'hWAT_D{day}_sem'])<10

    return hwat
