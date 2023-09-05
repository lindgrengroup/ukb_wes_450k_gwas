# +
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from adjustText import adjust_text

from ukb_wes_450k_gwas.utils import read_saige_gwas, obesity_phenotype_to_label_dict, LOCAL_PLOTS_DIR, get_effect_size_of_minor_allele

COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']

# def old1_get_saige_path(pheno, chrom, sex='both_sexes', assoc='variant'):
#     if assoc=='variant':
#         return f'{DATA_DIR}/old-02_saige_variant_test/saige_variant_test.{pheno}{f"-{sex}" if sex != "both_sexes" else ""}.chr{chrom}.tsv'

# def old2_get_saige_path(pheno, chrom, sex='both_sexes', assoc='gene_assoc'):
#     """
#     Run with variant and gene-based assocation simultaneously. Results are weird.
#     """
#     return f'{DATA_DIR}/02_saige_all_test/{sex}/saige_all_test.{pheno}{f"-{sex}" if sex != "both_sexes" else ""}.chr{chrom}.{assoc}_assoc.tsv'


def get_is_coding_nonsynonymous(df):
    return (
        df.consequence_category.notna()
        & ~df.consequence_category.isin({'non_coding', 'synonymous'})
    )

def clump_loci(df, clump_by='gene_symbol', sort_by='nlog10pval', ascending=False, verbose=False):
    df = df.sort_values(by=sort_by, ascending=ascending)
    
    n_variants = df.shape[0]
    loci_df_list = []
    
    if clump_by=='gene_symbol':
        max_iter = df.shape[0]
        i = 0
        while len(df)>0:
            top_hit = df.iloc[:1]

            chrom = top_hit['chrom'].values[0]
            position = top_hit['position'].values[0]
            gene_symbol = top_hit['gene_symbol'].values[0]
            
            is_in_locus = (
                (df['gene_symbol'] == gene_symbol)
                | (
                    (df['chrom'] == chrom)
                    & (df['position'] == position)
                  )
            )
            
            if verbose:
                print(f'chr{chrom}:{position} {gene_symbol} ({is_in_locus.sum()} variants)')
            
            top_hit['n_coding_nonsynonymous_in_locus'] = get_is_coding_nonsynonymous(df[is_in_locus]).sum()
            top_hit['n_variants_in_locus'] = is_in_locus.sum()
            
            df = df[~is_in_locus]
            
            loci_df_list.append(top_hit)
            i += 1
            
            if i > max_iter:
                raise RuntimeError('Max iterations reached during locus clumping')
            
    n_loci = len(loci_df_list)
    
    print(f'{n_variants} variants clumped into {n_loci} loci')
    
    return pd.concat(loci_df_list, axis=0)


def add_labels(x, y, labels, ax, highlight_labeled_points_kwargs={}, sig=False):
    labels_to_add = [
        ax.text(
            x[i], 
            y[i], 
            labels[i], 
        ) for i in range(len(x))
    ]

    if len(highlight_labeled_points_kwargs)>0:
        ax.scatter(
            x,
            y, 
            **highlight_labeled_points_kwargs
        )

    return labels_to_add

def plot_qq(df, nlog10pval_field = 'nlog10pval', maf_bins=None, logscale=False, 
    figsize=(9*1.2,6*1.2), dpi=100, title='', ax=None, scatter_kwargs={}, legend=True, print_maf_bin_counts=True,
    n_top_genes_to_label=0, clump_variants_into_loci=False, assoc='variant'):
    r'''Makes QQ plot for fields in `fields`.

    maf_bins should be a list of tuples, indicating min and max MAF per bin.
    Plots y-axis/x-axis in log scale if `logscale`=True.
    '''
    def get_expected(n):
        """Get expected -log10(p) for `n` observations
        """
        exp = -np.log10(np.linspace(start=1,stop=1/n,num=n)) # to account for weird kink at lower end in log-log scale: -np.log10(np.linspace(start=1-1/df_tmp.shape[0],stop=1/df_tmp.shape[0],num=df_tmp.shape[0]))

        return exp

    # Make 95% confidence interval
    def get_confidence_intervals(n, CI=0.95):
        k = np.arange(1,n+1)
        a = k
        b = n+1-k
        intervals=stats.beta.interval(CI, a, b) # get 95% CI with parameters a and b
        return intervals

    def get_lambda_gc(chisq_vec):
        return np.median(chisq_vec)/stats.chi2.ppf(q=0.5, df=1)

    df = df[df[nlog10pval_field].notna()].copy()
    
    if clump_variants_into_loci:
        raise ValueError('Not implemented yet')
        print('WARNING: Clumping into loci may artificially inflate the QQ plot because we use p-values from the lead (most-significant) variants to represent a locus')
        df = clump_loci(df)
        is_coding_nonsynonymous = df['n_coding_nonsynonymous_in_locus']>0
    elif assoc=='variant':
        is_coding_nonsynonymous = get_is_coding_nonsynonymous(df)

    n = df.shape[0]
    exp = get_expected(n)
    intervals = get_confidence_intervals(n)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.plot(exp,exp,'k--', label='Null')
    ax.fill_between(
        x=exp[::-1], # need to reverse order of array exp to match intervals
        y1=-np.log10(intervals[0]),
        y2=-np.log10(intervals[1]),
        color='lightgrey',
        label='95% CI'
    )

    lambda_gc_dict = {}

    if maf_bins is None:
        exp = get_expected(n=n)
        rank = np.sort(df[nlog10pval_field])
        ax.scatter(
            exp, 
            rank, 
            label='Observed', 
            **scatter_kwargs
        )
    else:
        if print_maf_bin_counts: print('MAF bin variant counts:')
        for i, (maf_min, maf_max) in enumerate(maf_bins):
            df_tmp = df[(df.maf>maf_min)&(df.maf<=maf_max)]
            maf_range_str = f'({maf_min}, {maf_max}]'
            if print_maf_bin_counts: 
                print(f'* {maf_range_str}: {len(df_tmp)} ({100*len(df_tmp)/len(df):.1f}%)')
            n_tmp = df_tmp[nlog10pval_field].notna().sum()
            exp_tmp = get_expected(n_tmp)
            rank_tmp = np.sort(df_tmp[nlog10pval_field])

            lambda_gc_tmp = get_lambda_gc(chisq_vec=df_tmp['chisq'])
            lambda_gc_dict[maf_range_str] = lambda_gc_tmp

            ax.scatter(
                exp_tmp,
                rank_tmp,
                color=plt.cm.viridis(i/(len(maf_bins)-0.5)),
                label=f'MAF:({maf_min}, {maf_max}]'+r', $\lambda_{GC}=$'+f'{lambda_gc_tmp:.2e}',
                **scatter_kwargs
            )

    if logscale:
        ax.set_yscale('log')
        ax.set_xscale('log')
        
        #    plt.ylim([min(-np.log10(intervals[0])),max(-np.log10(intervals[1]))])
    ax.set_xlabel(r'Expected -$\log_{10}(p)$')
    ax.set_ylabel(r'Observed -$\log_{10}(p)$')
    if legend:
        ax.legend(loc='upper left')
    
    gene_labels=[]
    if maf_bins is None and n_top_genes_to_label>0:
        gene_labels += add_labels(
            x=exp[-n_top_genes_to_label:],
            y=rank[-n_top_genes_to_label:],
            labels= df.sort_values(by=nlog10pval_field)['gene_symbol'].values[-n_top_genes_to_label:], 
            ax=ax, 
        )
    
    if len(gene_labels)>0:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim([ymin, ymax*1.2])
        adjust_text(
            texts=gene_labels, 
            ax=ax, 
            autoalign=True,
            arrowprops=dict(arrowstyle='->', color='black')
        )

    
    if assoc=='variant':
        lambda_gc_dict['all'] = get_lambda_gc(chisq_vec=df['chisq'])
        title += r'$\lambda_{GC}='+f'{lambda_gc_dict["all"]:.2e}$'
    ax.set_title(title)

    return df, lambda_gc_dict

def plot_manhattan(df, log_yscale=False, title='', ax=None, scatter_kwargs={}, nlog10pval_field='nlog10pval',
    sig_thresholds_to_plot=[], highlight_coding_nonsynonymous=False, colors=COLORS):
    
    start_bp = 0
    mid = []
    chrom_list = []
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    if 'chrom' in df.columns:
        chrom_field = 'chrom'
    else:
        chrom_field = 'CHR'
    for i, chrom in enumerate(range(1, 24)):
        max_pos=0
        
        if chrom==23: chrom='X'
#         try:
        gwas = df[df[chrom_field].astype(str)==str(chrom)]
    
            #         print(f'Variants with defined p-values chr{chrom}: {gwas.nlog10pval.notna().sum()}')
        try:
            if highlight_coding_nonsynonymous:
                tmp_df_dict = {}
                is_coding_nonsynonymous = gwas.consequence_category.isin({'other_missense','damaging_missense','pLoF'})
                tmp_df_dict['other'] = gwas[~is_coding_nonsynonymous]
                tmp_df_dict['coding_nonsynonymous'] = gwas[is_coding_nonsynonymous]
                for label in ['other', 'coding_nonsynonymous']:
                    ax.scatter(
                        tmp_df_dict[label]['position']+start_bp,
                        tmp_df_dict[label][nlog10pval_field], 
                        c=colors[i % 2] if label=='coding_nonsynonymous' else 'grey',
                        alpha=0.6 if label=='coding_nonsynonymous' else 0.1,
                        **scatter_kwargs
                    )
                    if len(tmp_df_dict[label]) > 0:
                        max_pos = max(max_pos, max(tmp_df_dict[label]['position']))
            else:
                ax.scatter(
                    gwas['position']+start_bp,
                    gwas[nlog10pval_field], 
                    c=colors[i % 2],
                    **scatter_kwargs
                )
                max_pos = max(max_pos, max(gwas['position']))
        except:
            max_pos = 0
            print(f'Failed chr{chrom}')
        mid.append(start_bp+(max_pos)/2)
        start_bp += max_pos
        start_bp += 5e7
        chrom_list.append(chrom)

    left, right = ax.get_xlim()
    for sig_threshold in sig_thresholds_to_plot:
        ax.plot([left, right], [-np.log10(sig_threshold)]*2, 'k--')
    ax.set_xticks(mid)
    ax.set_xticklabels(chrom_list)
    ax.set_title(title)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p)')
    ax.set_xlim([left, right])
    _, ymax = ax.get_ylim()
    if log_yscale:
        ax.set_yscale('symlog', linthresh=5)
        ax.set_ylim([0, ymax**1.05])

def plot_multiple_manhattan(df_dict, phenotype_group, pop, sex, pheno_list, ncols = 2, nrows=5):
    fig, axs = plt.subplots(figsize=(6*1.2*ncols,4*1.2*nrows), dpi=300, nrows=nrows, ncols=ncols)

    for i, pheno in enumerate(pheno_list):
        ax = axs[i//2][i%2]
        gwas_id = f'{pheno}-{sex}'
#         try:
        print(gwas_id)
        if gwas_id not in df_dict:
            df_dict[gwas_id] = read_saige_gwas(df_dict=df_dict, pheno=pheno, phenotype_group=phenotype_group, pop=pop, sex=sex)
        df = df_dict[gwas_id]
        plot_manhattan(
            df=df, 
            ax=ax,
            title=f'{obesity_phenotype_to_label_dict[pheno]}, sex={sex}\n'+f'(n={int(df.N.mean())})',
            scatter_kwargs={'s':5},
            highlight_coding_nonsynonymous=True
        )
        n_non_null = (df['p.value'].notna()).sum()
        print(f'* Variants with non-null p-values: {n_non_null}')
        print(f'* Variants passing GWS (pval<5e-8): {(df.nlog10pval > -np.log10(5e-8)).sum()}')
        # Nominal Bonferroni-corrected significance: (pval<0.05/[count of variants with non-null p-values])
#         nominal_bonf_thresh = 0.05/n_non_null
#         print(f'* Variants passing nominal Bonf.-corrected sig. (pval<{nominal_bonf_thresh:.3e}): {(df_dict[pheno].nlog10pval > -np.log10(nominal_bonf_thresh)).sum()}')
        print()
#         except:
#             print(f'Failed {pheno}')
    plt.tight_layout()

    fname=f'manhattan_plot.{"-".join(sorted(pheno_list))}.png'
    plt.savefig(f'{LOCAL_PLOTS_DIR}/{fname}', dpi=300)


def plot_maf_vs_effect(df, log_yscale=False, title='', ax=None, scatter_kwargs={}, 
    n_most_sig_variants_to_label_gene=0, n_highest_effect_variants_to_label_gene=0,clump_variants_into_loci=False,
    label_coding_nonsynonymous_only=True):
    
    df['abs_beta'] = df['BETA'].abs()
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    ax.set_title(title)
    ax.set_xlabel('MAF')
    ax.set_ylabel('abs(beta)')
    
    if clump_variants_into_loci:
        df = clump_loci(df)
    
    tmp_df_dict = {}
    if clump_variants_into_loci:
        is_coding_nonsynonymous = df['n_coding_nonsynonymous_in_locus']>0
    else:
        is_coding_nonsynonymous = df.consequence_category.isin({'other_missense','damaging_missense','pLoF'})
    tmp_df_dict['other'] = df[~is_coding_nonsynonymous]
    tmp_df_dict['coding_nonsynonymous'] = df[is_coding_nonsynonymous]
    
    point_size_scaling = 10
    
    for label in ['other', 'coding_nonsynonymous']:
        ax.scatter(
            tmp_df_dict[label]['maf'],
            tmp_df_dict[label]['abs_beta'], 
            c='tab:blue' if label=='coding_nonsynonymous' else 'grey',
            alpha=0.6 if label=='coding_nonsynonymous' else 0.1,
            s = point_size_scaling*tmp_df_dict[label]['nlog10pval'],
#             cmap='YlGnBu',
            **scatter_kwargs
        )
    
    gene_labels = []
    
    df_for_labels = tmp_df_dict['coding_nonsynonymous'] if label_coding_nonsynonymous_only else df
    
    def _add_gene_labels(sort_by, n_to_label):
        top_hits_df = df_for_labels.sort_values(by=sort_by, ascending=False).iloc[:n_to_label]
        
        highlight_labeled_points_kwargs = dict(
            facecolors='none',
            edgecolors='black',
            alpha=0.8,
            s = point_size_scaling*top_hits_df['nlog10pval'],
        )
        
        return add_labels(
            top_hits_df['maf'].values,
            top_hits_df['abs_beta'].values, 
            top_hits_df['gene_symbol'].values, 
            ax=ax, 
            highlight_labeled_points_kwargs=highlight_labeled_points_kwargs, 
        )
        
    if n_most_sig_variants_to_label_gene>0:
        gene_labels += _add_gene_labels(sort_by='nlog10pval', n_to_label=n_most_sig_variants_to_label_gene)
    if n_highest_effect_variants_to_label_gene>0:
        gene_labels += _add_gene_labels(sort_by='abs_beta', n_to_label=n_highest_effect_variants_to_label_gene)

    if len(gene_labels)>0:
        adjust_text(
            texts=gene_labels, 
            ax=ax, 
            autoalign=True,
            arrowprops=dict(arrowstyle='->', color='black')
        )


# +
def get_genes_with_pos(csq):
    csq = csq.copy()
    if 'gene_id' not in csq.columns:
        csq = csq.rename(columns={'Region':'gene_id'})
    csq['position'] = csq['MarkerID'].str.split(':', expand=True)[1].astype(int)
    genes_with_pos = csq[['gene_id','position']].groupby('gene_id').mean().reset_index()
    print(f'WARNING: This uses GRCh38 gene positions. Only use for GRCh37 (e.g. imputed v3 data) as an approximation')
    return genes_with_pos

def add_gene_position(df, genes_with_pos):
    if 'position' not in df.columns:
        if 'gene_id' not in genes_with_pos.columns:
            genes_with_pos = genes_with_pos.rename(columns={'Region':'gene_id'})
        df = df.merge(genes_with_pos[['gene_id','position']], on=['gene_id'], how='left')

    return df


# -

def plot_gene_manhattan(df, nlog10pval_field = 'nlog10Pvalue', min_maf=1e-6, log_yscale=False, title='',
    ax=None, n_top_genes_to_label=0, chrom_padding=5e7, colors=COLORS):
    
    start_bp = 0
    mid = []
    chrom_list = []
    plot_df_list = []

    if ax is None:
        print('creating subplots')
        fig, ax = plt.subplots(dpi=300)
        
    for i, chrom in enumerate(range(1, 24)):
        if chrom==23: chrom='X'

        gwas = df[df.CHR.astype(str)==str(chrom)]
        if len(gwas)==0:
            continue
        gwas['plot_position'] = gwas['position']+start_bp
    #         print(f'Missing gene position: {gwas.position.isna().sum()}')

        ax.scatter(
            gwas['plot_position'],
            gwas[nlog10pval_field], 
            c=colors[i % 2]
        )
        plot_df_list.append(gwas)
        mid.append(start_bp+(max(gwas['position']))/2)
        start_bp += max(gwas['position'])
        start_bp += chrom_padding # Space between chromosomes
        chrom_list.append(chrom)

    left, right = ax.get_xlim()
    ax.plot([left, right], [-np.log10(0.05/20e3)]*2, 'k--') # 20k genes
    #     plt.plot([left, right], [-np.log10(0.05/df['p.value'].notna().sum())]*2, 'k--', alpha=0.2)
    ax.set_xticks(mid)
    ax.set_xticklabels(chrom_list)
    ax.set_title(title)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p)')
    ax.set_xlim([left, right])
    _, ymax = ax.set_ylim()

    # Add labels for most significant genes
    plot_df = pd.concat(plot_df_list, axis=0)
    plot_df = plot_df.sort_values(by=nlog10pval_field, ascending=False)
    gene_labels=[]
    if n_top_genes_to_label>0:
        gene_labels += add_labels(
            x=plot_df['plot_position'].values[:n_top_genes_to_label],
            y=plot_df[nlog10pval_field].values[:n_top_genes_to_label],
            labels= plot_df.sort_values(by=nlog10pval_field, ascending=False)['gene_symbol'].values[:n_top_genes_to_label], 
            ax=ax, 
        )

    if len(gene_labels)>0:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim([ymin, ymax*1.2])
        adjust_text(
            texts=gene_labels, 
            ax=ax, 
            autoalign=True,
            arrowprops=dict(arrowstyle='->', color='black')
        )

    if log_yscale:
        ax.set_yscale('symlog', linthresh=5)
        ax.set_ylim([0, ymax**1.05])

    return df


def plot_miami(df1, df2, log_yscale=False, title='', scatter_kwargs={}):
    
    start_bp = 0
    mid = []
    chrom_list = []
    plt.figure(figsize=(12, 8))
    for chrom_i, chrom in enumerate(range(1, 24)):
        if chrom==23: chrom='X'
            
        max_pos=0
        for i, df in enumerate([df1, df2]):
            try:
                gwas = df[df.CHR==chrom]

    #             print(f'GWAS {i+1}: Variants with defined p-values chr{chrom}: {gwas.nlog10pval.notna().sum()}')

                plt.scatter(
                    gwas['position']+start_bp,
                     ((-1)**i)*gwas['nlog10pval'],  # flip second GWAS to be below
                    c=COLORS[chrom_i % 2],
                    **scatter_kwargs
                )
                max_pos = max(max_pos, max(gwas['position']))
            except:
                print(f'Failed: GWAS {i+1} chr{chrom}')
        mid.append(start_bp+(max_pos)/2)
        start_bp += max_pos
        start_bp += 5e7
        chrom_list.append(chrom)

    left, right = plt.xlim()
    for i in range(2):
        plt.plot([left, right], [-np.log10(5e-7)*(-1)**i]*2, 'k--')
        plt.plot([left, right], [-np.log10(0.05/df['p.value'].notna().sum())*(-1)**i]*2, 'k--', alpha=0.2)
    plt.xticks(mid, chrom_list)
    plt.title(title)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p)')
    plt.xlim([left, right])
    _, ymax = plt.ylim()
    if log_yscale:
        plt.yscale('symlog', linthresh=5)
        plt.ylim([0, ymax**1.05])

def plot_identity(ax):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.plot(*[[min([xmin, ymin]), min(xmax, ymax)]]*2, 'k--')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


def plot_field_comparison(df1, df2, field='nlog10pval', xlabel='Pheno1', ylabel='Pheno2', title=None,
    merge_on=['chrom','position','MarkerID','Allele1','Allele2'], suffixes=('_1','_2'),
    ax=None, scatter_kwargs={}, include_identity=True):
    
    both_phenos = df1.merge(df2, on=merge_on, suffixes=suffixes, how='outer')
    
    if '_signed' in field:
        for suffix in suffixes:
            sign = both_phenos['BETA'+suffix]/both_phenos['BETA'+suffix].abs()
            both_phenos[field+suffix] = (sign)*both_phenos[field.replace('_signed','')+suffix]
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12,8), dpi=300)
    ax.scatter(
        both_phenos[field+suffixes[0]],
        both_phenos[field+suffixes[1]],
        **scatter_kwargs
    )
    if include_identity:
        plot_identity(ax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title((field if title is None else title))
    
    return both_phenos

def create_custom_legend(pip_thresh, wes_pval_thresh, wes_maf_thresh, group_test, gene_pval_thresh):
    from matplotlib.lines import Line2D
    
    imputed = Line2D([0], [0], linestyle='', marker='o', label=f'Fine-mapped imputed data (PIP>={pip_thresh})', color='tab:blue')
    wes_variant = Line2D([0], [0], linestyle='', marker='o', label=f'WES coding non-syn ($P$<{wes_pval_thresh}, MAF<{wes_maf_thresh})', color='tab:orange')
    wes_gene = Line2D([0], [0], linestyle='', marker='o', label=f'Gene burden test ({group_test} $P$<{gene_pval_thresh})', color='tab:green')

    plt.legend(handles=[imputed, wes_variant, wes_gene])
    
def plot_trumpet(df_dict, gwas_id, pip_thresh = 0.99, wes_pval_thresh = 5e-8, wes_maf_thresh = 0.01, 
                 gene_pval_thresh = 0.05/20e3, group_test='pLoF;damaging_missense', pval_burden_field='Pvalue_Burden', 
                 beta_burden_field='BETA_Burden', nlog10pval_burden_field='nlog10Pvalue_Burden',
                include_legend=True, figsize=(6,6), title=''):
    imputed_gwas_id = f'{gwas_id}-imputed'
    imputed_gwas = df_dict[imputed_gwas_id]
    imputed_gwas['BETA_minor_allele'] = get_effect_size_of_minor_allele(imputed_gwas)
    
    wes_gwas_id = f'{gwas_id}-wes'
    wes_gwas = df_dict[wes_gwas_id]
    wes_gwas['BETA_minor_allele'] = get_effect_size_of_minor_allele(wes_gwas)
    
    imp = imputed_gwas[imputed_gwas.prob>=pip_thresh]

    wes = wes_gwas[
        (wes_gwas['p.value']<wes_pval_thresh)
        & (wes_gwas['maf']<wes_maf_thresh)
        & is_coding_nonsynonymous(wes_gwas)
    ]
    wes['abs_BETA'] = wes['BETA']
    # In each gene, only keep highest effect size significant variant 
    wes = wes.sort_values(by='abs_BETA', ascending=False).drop_duplicates(subset='gene_symbol', keep='first')
    
    max_maf = 0.01
    df_dict, genes = get_gene_aggregated_maf(df_dict=df_dict, gwas_id=gwas_id, group_test=group_test, max_maf=max_maf)
    burden = genes[
        (genes['Group']==group_test)
        &(genes['max_MAF']==max_maf)
        &(genes[pval_burden_field]<gene_pval_thresh)
    ]

    fig, ax = plt.subplots(figsize=figsize, dpi=300)

    for df in [imp, wes]:
        ax.scatter(
            df['maf'],
            df['BETA_minor_allele'],
            s=3*df['nlog10pval'],
            alpha=0.8,
        )
        
    ax.scatter(
        burden['MAF_aggregated'],
        burden[beta_burden_field],
        s=3*burden[nlog10pval_burden_field],
        alpha=0.8,
    )

    ax.set_xscale('log')
    ax.set_xlabel('MAF')
    ax.set_ylabel('Effect size of minor allele\n(burden effect size for gene-level results)')
    ax.set_title(title)
    if include_legend:
        create_custom_legend(
            pip_thresh=pip_thresh,
            wes_pval_thresh=wes_pval_thresh,
            wes_maf_thresh=wes_maf_thresh,
            group_test=group_test,
            gene_pval_thresh=gene_pval_thresh
        )

    gene_labels = []

    gene_labels += add_labels(
        wes['maf'].values,
        wes['BETA_minor_allele'].values, 
        wes['gene_symbol'].values, 
        ax=ax
    )
    
    gene_labels += add_labels(
        burden['MAF_aggregated'].values,
        burden[beta_burden_field].values, 
        burden['gene_symbol'].values, 
        ax=ax
    )

    adjust_text(
        texts=gene_labels, 
        ax=ax, 
        autoalign=True,
        arrowprops=dict(arrowstyle='-', color='black', alpha=0.5)
    )
    
    plt.show()
    
    
def plot_regeneron_trumpet(df_dict, phenotype_group, gwas_id, pip_thresh = 0.90, imputed_data_pval_thresh=5e-8, 
                 gene_pval_thresh = 0.05/20e3, group_tests=['pLoF','pLoF;damaging_missense'], pval_burden_field='Pvalue_Burden', 
                 max_maf = 0.01, beta_burden_field='BETA_burden_duncan', nlog10pval_burden_field='nlog10Pvalue_Burden',
                 include_legend=True, figsize=(6,6), title='', include_gene_labels=True, yticks_list=[], n_imputed_gene_labels=0,
                 beta_scaling=1, ylabel='Effect size', axis_zoom=0, dpi=500, savefig=False):
    
    imputed_gwas_id = f'{gwas_id}-imputed'
    imputed_gwas = df_dict[imputed_gwas_id]
    imputed_gwas['beta_minor_allele'] = get_effect_size_of_minor_allele(imputed_gwas)
    
    imp = imputed_gwas[
        (imputed_gwas['prob']>=pip_thresh)
        & (imputed_gwas['pval']<imputed_data_pval_thresh)
    ]
    
    if gwas_id.split('-')[0]=='body_mass_index_bmi':
        # Use new version of phenotype name
        gwas_id = gwas_id.replace('body_mass_index_bmi', 'bmi')
        
    df_dict, genes = get_gene_aggregated_maf_and_burden_beta(
        df_dict=df_dict, 
        phenotype_group=phenotype_group,
        gwas_id=gwas_id, 
        max_maf=max_maf, 
        csq_categories_list=[group.split(';') for group in group_tests]
    )
#     print(genes.columns)
    burden = genes[
        (genes['Group'].isin(group_tests))
        &(genes['max_MAF']==max_maf)
        &(genes[pval_burden_field]<gene_pval_thresh)
    ]
    
    # For genes with multiple significant results, keep the result with highest effect size
    burden = burden.sort_values(by=beta_burden_field, ascending=False)
    burden = burden.drop_duplicates(subset='gene_symbol')
    
#     burden['gene_symbol'] = burden['gene_symbol']+'-'+burden['Group'].str.replace('damaging_missense','dm')

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_axisbelow(True)
    ax.grid(c='lightgrey', alpha=0.3)
    
    # Plot imputed results
    ax.scatter(
        imp['maf'],
        beta_scaling*imp['beta_minor_allele'],
        s=30+800*imp['beta_minor_allele'].abs(),
#         alpha=0.8,
        c='#4f9db9'
    )

    # Plot burden test
    ax.scatter(
        burden['burden_maf'],
        beta_scaling*burden[beta_burden_field],
        s=30+800*burden[beta_burden_field].abs(), #3*burden[nlog10pval_burden_field],
#         alpha=0.8,
        c='#f3641f'
    )

    fontname=None#'Helvetica Neue'
    fontweight=None #100
    ax.set_xscale('log')
    ax.set_xlabel('Minor allele frequency', fontname=fontname, fontweight=fontweight) #fontweight=600
    ax.set_ylabel(ylabel, fontname=fontname, fontweight=fontweight)
    ax.tick_params(labelsize=8)
    if len(yticks_list)>0:
        ax.set_yticks(ticks=yticks_list)
        ax.set_yticklabels(labels=yticks_list)
    ax.set_title(title, fontname=fontname, fontweight=fontweight)
    if include_legend:
        create_custom_legend(
            pip_thresh=pip_thresh,
            wes_pval_thresh=wes_pval_thresh,
            wes_maf_thresh=wes_maf_thresh,
            group_test='+'.join(group_tests),
            gene_pval_thresh=gene_pval_thresh
        )

    if include_gene_labels:
        gene_labels = []

        gene_labels += add_labels(
            burden['burden_maf'].values,
            beta_scaling*burden[beta_burden_field].values, 
            burden['gene_symbol'].values, 
            ax=ax
        )
        
        if n_imputed_gene_labels>0:
            imp['abs_beta'] = imp['beta'].abs()
            imp = imp.sort_values(by='abs_beta', ascending=False)
            gene_labels += add_labels(
                imp['maf'].values[:n_imputed_gene_labels],
                beta_scaling*imp['beta_minor_allele'].values[:n_imputed_gene_labels], 
                imp['gene_symbol'].str.split('/', expand=True)[0].values[:n_imputed_gene_labels], 
                ax=ax
            )

        adjust_text(
            texts=gene_labels, 
            ax=ax, 
            fontstyle='italic',
            autoalign=True,
            arrowprops=dict(arrowstyle='-', color='black', alpha=0.5)
        )
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ylim = plt.ylim()
    ax.set_ylim([y*1.05 for y in ylim])
    xlim = plt.xlim()
    ax.set_xlim([0.8*xlim[0], xlim[1]])
    plt.tight_layout()
    
    if savefig:    
        pass
#         plt.savefig(f'/Users/nbaya/Downloads/trumpet.{gwas_id}.pdf', dpi=dpi)
    
    plt.show()
    
    
if __name__=='__main__':
    pass
