# +
import pandas as pd

GWS_THRESHOLD = 5e-8 # Genome-wide significance threshold
WINDOW_WIDTH_BP = int(3e6) # 3 Mb


# -

def get_lead_variants(sig_loci, sort_by='nlog10pval_outlier', ascending=False, verbose=True):
    """
    Get most significant variants, excluding variants within a window around the lead variant
    """
    MHC_START=25e6
    MHC_STOP=34e6

    in_mhc = (
            (str(sig_loci.chrom) == '6')
            & (sig_loci.position >= MHC_START)
            & (sig_loci.position <= MHC_STOP)
        )

    if in_mhc.sum() > 0:
        sig_loci = sig_loci.loc[~in_mhc,:]
        if verbose: print(f'* Excluding {in_mhc.sum()} variants in MHC (chr6:{MHC_START}-{MHC_STOP})')

    sig_loci = sig_loci.sort_values(by=sort_by, ascending=ascending) # Highest value (most significant if using -log10(pval)) first

    lead_variant_row_list = []
    
    while len(sig_loci)>0:
        lead_variant_row = sig_loci.iloc[:1,:]
        
        lead_variant_row_list.append(lead_variant_row)
        
        lead_variant = lead_variant_row.iloc[0,:]

        if 'gene_id' in sig_loci.columns:
            lead_variant_gene = lead_variant_row['gene_id'].values[0]
            # NOTE: 
            is_in_same_gene = sig_loci['gene_id'] == lead_variant_gene
        else:
            is_in_same_gene = True
        
        window_start = lead_variant.position - WINDOW_WIDTH_BP/2
        window_stop = lead_variant.position + WINDOW_WIDTH_BP/2
        
        in_window = (
            (sig_loci.chrom==lead_variant.chrom)
            & (sig_loci.position>=window_start)
            & (sig_loci.position<=window_stop)
        ) | is_in_same_gene

        if verbose: print(f'* Variants in window around chr{lead_variant.chrom}:{lead_variant.position}:{lead_variant.Allele1}:{lead_variant.Allele2}: {in_window.sum()}')
        if 'gene_id' in sig_loci.columns:
            if verbose: print(f'  - Variants in same gene: {is_in_same_gene.sum()}')
        
        sig_loci = sig_loci.loc[~in_window,:]
        
    lead_variants = pd.concat(lead_variant_row_list)
    lead_variants['window_start'] = lead_variants.position-int(WINDOW_WIDTH_BP/2)
    lead_variants['window_stop'] = lead_variants.position+int(WINDOW_WIDTH_BP/2)
    lead_variants = lead_variants.sort_values(by=['chrom', 'position'])
    
    return lead_variants

def merge_overlapping_windows(lead_variants, variant_id_field='MarkerID'):
    merged_lead_variants_list = []

    # WARNING: The following code requires that variants are sorted by genetic position
    lead_variants = lead_variants.sort_values(by=['chrom', 'position'])

    while len(lead_variants)>0:
        # Get first row from table of lead variants
        variant_row = lead_variants.iloc[:1,:]
        variant = variant_row.iloc[0,:]
        
        merged_lead_variants_list.append(variant_row)

        # Remove first row
        lead_variants = lead_variants.iloc[1:,:]

        # Get lead variants on the same chrom that are different from the current lead variant
        same_chrom = lead_variants[
            (lead_variants.chrom==variant.chrom)
            & (lead_variants[variant_id_field]!=variant[variant_id_field])
        ]

        if len(same_chrom)>0:
            merged_variant_ids = []

            for i, other_variant in same_chrom.iterrows():

                # Boolean indicating whether other variant has window overlapping with variant
                is_overlapping = (
                    ( # Other variant window to the right and overlapping
                        (other_variant.window_start>variant.window_start)
                        & (other_variant.window_start<=variant.window_stop)
                    ) 
                    | ( # other variant window to the left and overlapping
                        (other_variant.window_start<variant.window_start)
                        & (other_variant.window_stop>=variant.window_start)
                    )
                )

                if is_overlapping:
                    # Redefine window start and stop positions
                    variant_row['window_start'] = min([variant.window_start, other_variant.window_start])
                    variant_row['window_stop'] = max([variant.window_stop, other_variant.window_stop])
                    merged_variant_ids.append(other_variant[variant_id_field])

            was_merged = lead_variants[variant_id_field].isin(merged_variant_ids)
            
            # Remove lead variants which have been merged into the current locus
            lead_variants = lead_variants[~was_merged]  

    merged_lead_variants = pd.concat(merged_lead_variants_list, axis=0)
    
    return merged_lead_variants
