import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']
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

# The four phenotypes with high sample size (>300 samples)
high_sample_size_obesity_phenos=[
    'body_mass_index_bmi',
    'whr_adj_bmi',
    'bmi_impedance',
    'body_fat_percentage'
]
# -

BIOMARKERS = pd.read_csv(
    '/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas/data/phenotypes/phenotype_list.qced_biomarkers.txt', 
    sep='\t', 
    names=['phenotype']
).phenotype.tolist()

SET_TEST_GROUPS = [
    'pLoF',
    'damaging_missense',
    'synonymous',
    'pLoF;damaging_missense',
    'pLoF;damaging_missense;synonymous',
    'Cauchy'   
]

# +
# pheno_to_label_dict = dict(zip(BIOMARKERS, BIOMARKERS)) # placeholder
# # pheno_to_label_dict = dict(zip(BIOMARKERS, [p.replace('_',' ') for p in BIOMARKERS]))
# # pheno_to_label_dict['c_reactive_protein'] = 'C-reactive protein'
# -

DNAX_DATA_DIR='/mnt/project/saige_pipeline/data'
LOCAL_WD='/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas'
LOCAL_DATA_DIR=f'{LOCAL_WD}/data'
LOCAL_PLOTS_DIR=f'{LOCAL_WD}/plots'

DATA_DIR=LOCAL_DATA_DIR

def get_saige_results_path(pheno, chrom, phenotype_group='qced_biomarkers', 
    pop='eur', sex='both_sexes', assoc='variant'):
    return f'{DATA_DIR}/02_saige_{assoc}_test/{phenotype_group}/{pop}/{sex}/saige_{assoc}_test.{pheno}-{pop}-{sex}.chr{chrom}.tsv.gz'

def get_consequence_path(chrom):
    return f'{DATA_DIR}/annotations/wes_450k/ukb_wes_450k.qced.chr{chrom}.worst_csq_by_gene_canonical.txt.gz'

def get_consequence_file(chrom):
    return pd.read_csv(get_consequence_path(chrom), compression='gzip', delim_whitespace=True)

def get_all_consequence_files():
    """Only autosomes
    """
    csq_df_list = [get_consequence_file(chrom=chrom) for chrom in range(1,23)]
    csq = pd.concat(csq_df_list, axis=0)
    csq = csq.rename(columns={'varid':'MarkerID'})
    return csq

def get_unique_variant_genes_and_consequences_path():
    return f'{DATA_DIR}/annotations/wes_450k/ukb_wes_450k.qced.worst_csq_by_gene_canonical.grouped_by_variant.tsv.gz'


def write_unique_variant_genes_and_consequences():
    def concat_genes(x):
        return '/'.join(x)
    
    csq = get_all_consequence_files()
    csq = csq.drop(columns=['locus','alleles'])
    csq_grouped = csq.groupby(['MarkerID','consequence_category']).agg((concat_genes))
    
    csq_grouped.to_csv(
        get_unique_variant_genes_and_consequences_path(),
        compression='gzip',
        sep='\t',
        index=False
    )


def get_unique_variant_genes_and_consequences():
    csq_grouped = pd.read_csv(
        get_unique_variant_genes_and_consequences_path(),
        sep='\t',
        compression='gzip'
    )
    
    return csq_grouped


def get_regenie_path(pheno, chrom, test='additive', outlier_type=None):
    pheno = pheno.replace('_','-')
    if outlier_type is None:
        return f'/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas/data/regenie/qced_biomarkers/eur/both_sexes/step2_{test}_chr{chrom}_{pheno}_qced.regenie.gz'

def read_saige_gwas(df_dict, pheno, phenotype_group, pop='eur', sex='both_sexes'):

    df_list = []

    for chrom in list(range(1,23))+['X']:
        try:
            path = get_saige_results_path(
                phenotype_group=phenotype_group,
                pheno=pheno, 
                sex=sex, 
                pop=pop,
                assoc='variant',
                chrom=chrom,
            )
            df_tmp = pd.read_csv(path, sep='\t', compression='gzip')
            df_list.append(df_tmp)
        except:
            print(path)
            print(f'No file for {pheno} pop={pop} sex={sex} chr{chrom}')
    df =  pd.concat(df_list, axis=0)
    df = df.rename(columns={
        'POS':'position',
        'CHR':'chrom'
    })
    df['CHISQ'] = (df.BETA/df.SE)**2
    df['p.value'] = df['p.value'].astype(float)
    df['nlog10pval'] = -1*np.log10(df['p.value'])
    df.loc[df['p.value']==0, 'nlog10pval'] = 315 # Use as placeholder for "very significant"
    
    df['maf'] = pd.concat([df['AF_Allele2'], 1-df['AF_Allele2']], axis=1).min(axis=1)
#     print(f'Removing {(df.CHR!="UR").sum()} lines with "UR"')
#     df = df[df.CHR!='UR']
    df['position'] = df['position'].astype(int)
    csq = df_dict['csq'][['MarkerID','gene_id','gene_symbol','consequence_category']]
    
    df = df.merge(csq, on='MarkerID', how='left')
    
    return df

def read_saige_gene_assoc(df_dict, pheno, phenotype_group, pop='eur', sex='both_sexes'):

    df_list = []

    for chrom in list(range(1,23))+['X']:
        path = get_saige_results_path(pheno=pheno, phenotype_group=phenotype_group, chrom=chrom, sex=sex, assoc='set')
        try:
            df_tmp = pd.read_csv(path, sep='\t')
            df_tmp['CHR'] = str(chrom)
            for pval_field in ['Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT']:
                df_tmp[f'nlog10{pval_field}'] = -np.log10(df_tmp[pval_field])
            df_list.append(df_tmp)
        except:
            print(f'No file for {pheno} sex={sex} chr{chrom}: {path}')
    df =  pd.concat(df_list, axis=0)
    
    gene_id_to_symbol_df = df_dict['csq'][['gene_symbol','gene_id']].drop_duplicates(subset='gene_id')
    df = df.merge(gene_id_to_symbol_df, left_on='Region', right_on='gene_id', how='left')
    
    return df


def read_regenie_results(pheno):
    get_path = lambda chrom: get_regenie_path(pheno=pheno, chrom=chrom)
    df_list = [pd.read_csv(get_path(chrom), compression='gzip', sep=' ') for chrom in range(1,23)]
    df = pd.concat(df_list, axis=0)

    df = df.rename(columns={
        'CHROM':'chrom',
        'GENPOS':'position',
        'ID':'MarkerID',
        'LOG10P':'nlog10pval',
        'ALLELE0':'Allele1',
        'ALLELE1':'Allele2',
        'A1FREQ':'AF_Allele2',
    })
    
    return df

def read_genebass_results(pheno, assoc='variant'):
    genebass_suffix = genebass_file_dict[pheno]
    if assoc=='gene':
        gb = pd.read_csv(f'{DATA_DIR}/genebass/gene-burden-results-exomes_pLoF_continuous-{genebass_suffix}.csv')
    elif assoc=='variant':
        gb = pd.read_csv(f'{DATA_DIR}/genebass/single-variant-associations-exomes_continuous-{genebass_suffix}.csv')
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
