# +
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from adjustText import adjust_text

from ukb_wes_450k_gwas.utils import DATA_DIR
# -
# def old1_get_saige_path(pheno, chrom, sex='both_sexes', assoc='variant'):
#     if assoc=='variant':
#         return f'{DATA_DIR}/old-02_saige_variant_test/saige_variant_test.{pheno}{f"-{sex}" if sex != "both_sexes" else ""}.chr{chrom}.tsv'

# def old2_get_saige_path(pheno, chrom, sex='both_sexes', assoc='gene_assoc'):
#     """
#     Run with variant and gene-based assocation simultaneously. Results are weird.
#     """
#     return f'{DATA_DIR}/02_saige_all_test/{sex}/saige_all_test.{pheno}{f"-{sex}" if sex != "both_sexes" else ""}.chr{chrom}.{assoc}_assoc.tsv'


# +
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


# -

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
        color='k',
        alpha=0.1,
        label='95% CI'
    )

    lambda_gc_dict = {}

    if maf_bins is None:
        exp = get_expected(n=n)
        rank = np.sort(df[nlog10pval_field])
        ax.scatter(exp, rank, label='Observed', **scatter_kwargs)
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

            lambda_gc_tmp = get_lambda_gc(chisq_vec=df_tmp['CHISQ'])
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
        lambda_gc_dict['all'] = get_lambda_gc(chisq_vec=df['CHISQ'])
        title += r'$\lambda_{GC}='+f'{lambda_gc_dict["all"]:.2e}$'
    ax.set_title(title)

    return df, lambda_gc_dict

def plot_manhattan(df, log_yscale=False, title='', ax=None, scatter_kwargs={}, 
    sig_thresholds_to_plot=[], highlight_coding_nonsynonymous=False):
    
    start_bp = 0
    mid = []
    chrom_list = []
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))
        
    for i, chrom in enumerate(range(1, 24)):
        max_pos=0
        
        if chrom==23: chrom='X'
#         try:
        gwas = df[df.CHR.astype(str)==str(chrom)]
    
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
                        tmp_df_dict[label]['nlog10pval'], 
                        c=COLORS[i % 2] if label=='coding_nonsynonymous' else 'grey',
                        alpha=0.6 if label=='coding_nonsynonymous' else 0.1,
                        **scatter_kwargs
                    )
                    if len(tmp_df_dict[label]) > 0:
                        max_pos = max(max_pos, max(tmp_df_dict[label]['position']))
            else:
                ax.scatter(
                    gwas['position']+start_bp,
                    gwas['nlog10pval'], 
                    c=COLORS[i % 2],
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


def plot_gene_manhattan(df, nlog10pval_field = 'nlog10Pvalue', min_maf=1e-6, log_yscale=False, title=''):
    df = df.merge(csq, on=['Region', 'CHR'], how='left')
    
    start_bp = 0
    mid = []
    chrom_list = []
    plt.figure(figsize=(12, 8))
    for i, chrom in enumerate(range(1, 23)):
        if chrom==23: chrom='X'
            
        gwas = df[df.CHR.astype(str)==str(chrom)]
#         print(f'Missing gene position: {gwas.position.isna().sum()}')
        
        plt.scatter(
            gwas['position']+start_bp,
            gwas[nlog10pval_field], 
            c=COLORS[i % 2]
        )
        mid.append(start_bp+(max(gwas['position']))/2)
        start_bp += max(gwas['position'])
        start_bp += 5e7
        chrom_list.append(chrom)

    left, right = plt.xlim()
    plt.plot([left, right], [-np.log10(0.05/20e3)]*2, 'k--') # 20k genes
#     plt.plot([left, right], [-np.log10(0.05/df['p.value'].notna().sum())]*2, 'k--', alpha=0.2)
    plt.xticks(mid, chrom_list)
    plt.title(title)
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p)')
    plt.xlim([left, right])
    _, ymax = plt.ylim()
    if log_yscale:
        plt.yscale('symlog', linthresh=5)
        plt.ylim([0, ymax**1.05])

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

def plot_field_comparison(df1, df2, field='nlog10pval', xlabel='Pheno1', ylabel='Pheno2', title=None,
    merge_on=['CHR','position','MarkerID','Allele1','Allele2'], suffixes=('_1','_2'),
    ax=None, scatter_kwargs={}, plot_identity=True):
    
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
    if plot_identity:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.plot(*[[min([xmin, ymin]), min(xmax, ymax)]]*2, 'k--')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title((field if title is None else title))
    
    return both_phenos


if __name__=='__main__':
    pass

