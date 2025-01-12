#!/usr/bin/env python3

import argparse
import copy
import os
import re
import warnings

from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, euclidean, cityblock
from scipy.cluster.hierarchy import linkage
import seaborn as sns

EPSILON = np.finfo(np.float64).resolution
CHR_ORDER = dict({f'{i}': i for i in range(1, 23, 1)}, **{'X': 23, 'Y': 24})
BASES = ['A', 'C', 'G', 'T']
VCF_COLS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
    'FORMAT', 'sample']
GT_MAP = {'0/0': 0, '0/1': 1, '1/1': 2, './.': np.nan}
PID2SID = {'PID1677': 'SID3841','PID1712': 'SID3867', 'PID2178': 'SID5962'}
SID2PAPER = {'SID3841': 'S1', 'SID3867': 'S2', 'SID5962': 'S3'}
Cluster2PAPER = {'0': 'MS1:S3', '1': 'MS1:S1', '2': 'MS1:S2'}

VAF_cmap = LinearSegmentedColormap.from_list('my_gradient', (
    (0.000, (1.000, 1.000, 1.000)),
    (0.500, (1.000, 0.616, 0.000)),
    (1.000, (1.000, 0.000, 0.000)))
)

FONTSIZE = 30
DPI = 300
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'lines.linewidth': 2,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'large',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 50
})
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['xtick.bottom'] = True


def get_rna_profiles(in_files, dna_loci):
    matched_map = {}
    for in_file in in_files:
        df_in = pd.read_csv(in_file, comment='#', sep='\t', header=None)
        df_in.columns = VCF_COLS
        df_in.rename({'CHROM': 'Chr', 'POS': 'Position'}, axis=1, inplace=True)
        df_in['Chr'] = df_in['Chr'].str.replace('chr', '')
        df_in.set_index(['Chr', 'Position'], inplace=True)
        overlap = df_in.index.isin(dna_loci)
        df_in = df_in.loc[overlap]
        df_in = df_in[df_in['FILTER'] != 'BACKGROUND']

        name = os.path.basename(in_file).split("_")[0]
        try:
            SID = PID2SID[name]
        except KeyError:
            name = f'{name}_RNA'
        else:
            try:
                name = f'{SID2PAPER[SID]}_RNA'
            except KeyError:
                name = f'{SID}_RNA'
        matched_map[name] = in_file

        # souporcell output: GT:AO:RO:T:E:GO:GN
        # INFO=<ID=RO,Description="Count of full observations of the reference haplotype.">
        # INFO=<ID=AO,Description="Count of full observations of this alternate haplotype.">
        alt = df_in['sample'].apply(lambda x: int(x.split(':')[1]))
        ref = df_in['sample'].apply(lambda x: int(x.split(':')[2]))
        dp = ref + alt
        VAF = np.clip((alt + EPSILON) / (dp + EPSILON), EPSILON, 1 - EPSILON)

        new_df = VAF.to_frame(name=name)
        new_df[f'snp_{name}'] = df_in['REF'] + '>' + df_in['ALT']

        try:
            if name in df:
                df_samples = df[[name]] \
                    .merge(new_df[name], left_index=True, right_index=True, how='outer')
                df[name] = df_samples.mean(axis=1)
            else:
                df = df.merge(new_df, left_index=True, right_index=True, how='outer')
        except NameError:
            df = new_df

    sample_cols = [i for i in df.columns if not i.startswith('snp_')]
    print(f'#SNPs:\n{(df[sample_cols] > 0.0).sum()}')
    return df, matched_map


def get_snp_array_profiles(in_file, dna_loci):
    matched_map = {}

    df_in = pd.read_csv(in_file, dtype={0: str, 3: str}, sep='\t', skiprows=9)
    df_in['Chr'] = df_in['Chr'].str.replace('chr', '')
    df_in.set_index(['Chr', 'Position'], inplace=True)
    overlap = df_in.index.isin(dna_loci)
    df_in = df_in.loc[overlap]
    # Remove not-standard basepair substitutions
    df_in = df_in[(df_in['Allele1 - Forward'].isin(BASES)) \
        & (df_in['Allele2 - Forward'].isin(BASES))]
    # Remove SNPs flagged as Ilumina Duplicate if 'normal' SNP is also present
    df_in = df_in[~(df_in['SNP Name'].str.contains('IlmnDup') \
        & df_in.index.duplicated(keep=False))]

    for sample, df_s in df_in.groupby('Sample Name'):
        try:
            if sample in SID2PAPER:
                SID = sample
            else:
                SID = PID2SID[sample]
        except KeyError:
            sample = f'{sample}_SNParray'
        else:
            try:
                sample = f'{SID2PAPER[SID]}_SNParray'
            except KeyError:
                sample = f'{SID}_SNParray'
        matched_map[sample] = in_file

        new_df = df_s['B Allele Freq'].to_frame(name=sample)
        new_df[f'snp_{sample}'] = df_s['Allele1 - Forward'] \
            + '>' + df_s['Allele2 - Forward']
        try:
            df_out = df_out.merge(new_df,
                left_index=True, right_index=True, how='outer')
        except NameError:
            df_out = new_df
    
    sample_cols = [i for i in df_out.columns if not i.startswith('snp_')]
    print(f'#SNPs:\n{(df_out[sample_cols] > 0.0).sum()}')
    return df_out, matched_map


def get_dna_profiles(in_files):
    samples = pd.DataFrame([], columns=['scDNA_file'])
    samples.index.name = 'scDNA_name'

    FILE_PATTERN = r'(\d+)(\.whitelist)?\.filtered_variants\.csv'
    for i, in_file in enumerate(in_files):
        if re.search(FILE_PATTERN, in_file):
            cluster = re.search(FILE_PATTERN, in_file).group(1)
            name = f'{Cluster2PAPER[cluster]}_DNA'
        else:
            name = f'{i}_DNA'

        samples.loc[name, 'scDNA_file'] = in_file

        df_new = get_dna_sample_profile(in_file, name)
        try:
            df = df.merge(df_new,
                left_index=True, right_index=True, how='outer')
        except NameError:
            df = df_new
    return df, samples


def get_dna_sample_profile(in_file, name):
    df = pd.read_csv(in_file, dtype={'CHR': str})
    snp_str = df['REF'].fillna('*') + '>' + df['ALT'].fillna('*')

    try:
        ref = df.iloc[:,7:].map(lambda x: int(x.split(':')[0]))
        alt = df.iloc[:,7:].map(lambda x: int(x.split(':')[1]))
    except AttributeError: # Older pandas versions < 2.1.0
        ref = df.iloc[:,7:].applymap(lambda x: int(x.split(':')[0]))
        alt = df.iloc[:,7:].applymap(lambda x: int(x.split(':')[1]))
    dp = ref + alt
    VAF = np.clip((alt + EPSILON) / (dp + EPSILON), EPSILON, 1 - EPSILON)

    df_out = pd.DataFrame(np.average(VAF, axis=1, weights=dp), columns=[name])
    df_out[f'snp_{name}'] = snp_str
    df_out[['Chr', 'Position']] = df[['CHR', 'POS']]
    df_out.set_index(['Chr', 'Position'], inplace=True)

    return df_out


def plot_profiles(data, out_file):
    dist_row = pdist(data.values, euclidean_nan)
    Z_row = linkage(np.nan_to_num(dist_row, 1), 'ward')

    cm = sns.clustermap(
        data,
        row_linkage=Z_row,
        row_cluster=True,
        col_cluster=False,
        col_colors=None,
        row_colors=None,
        vmin=0, center=0.5, vmax=1,
        cmap=VAF_cmap,
        # dendrogram_ratio=(0.1, 0.1),
        figsize=(25, 15),
        cbar_kws={'shrink': 0.5, 'ticks': [0, 0.5, 1]},
        cbar_pos=(0.15, 0.8, 0.01, 0.075),
        tree_kws=dict(linewidths=2),
        xticklabels=True
    )
    hm = cm.ax_heatmap
    hm.set_facecolor('#636161')

    hm.set_ylabel('\nSamples', fontsize=FONTSIZE)
    hm.set_yticklabels(cm.ax_heatmap.get_yticklabels(), rotation='horizontal')

    hm.set_ylabel('SNPs', fontsize=FONTSIZE)
    hm.set_xticklabels(hm.get_xticklabels(), fontsize=FONTSIZE/1.5, va='top')

    cm.ax_cbar.set_title('VAF', fontsize=FONTSIZE)

    cm.fig.tight_layout()
    print(f'Saving SNP profiles heatmap: {out_file}')
    cm.fig.savefig(out_file, dpi=DPI)


def plot_distance(df_in, out_file):
    df = copy.deepcopy(df_in)

    col_no = 1
    row_no = 1
    fig, ax = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no*12, row_no*12))

    sns.set(font_scale=2.5)

    df.sort_index(ascending=False, inplace=True)
    df = df[sorted(df.columns)]

    df.columns = [i.split('_')[0] for i in df.columns]
    df.index = [i.split('_')[0] for i in df.index]

    hm = sns.heatmap(
        df,
        annot=True,
        square=True,
        cmap='viridis_r',
        cbar_kws={'shrink': 0.5, 'label': 'Distance'},
        linewidths=0,
        ax=ax
    )

    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=FONTSIZE)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=FONTSIZE)
    ax.set_ylabel(df_in.index.name, fontsize=FONTSIZE+10)
    ax.set_xlabel(df_in.columns.name, fontsize=FONTSIZE + 10)

    fig.tight_layout()
    print(f'Saving sample distance heatmap: {out_file}')
    fig.savefig(out_file, dpi=300)


def assign_clusters(df_dist, sample_df, matched_dict, assignment_output):
    df = copy.deepcopy(df_dist)
    for i in range(df_dist.shape[0]):
        cl_idx, other_idx = np.where(df == df.min().min())
        cl = df.index[cl_idx[0]]
        other = df.columns[other_idx[0]]

        print(f'Assigning: {cl} -> {other}')
        sample_df.loc[cl, 'matched_name'] = other
        sample_df.loc[cl, 'matched_file'] = matched_dict[other]
        
        df.drop(cl, axis=0, inplace=True)
        df.drop(other, axis=1, inplace=True)

    print(f'Saving sample assignments: {assignment_output}')
    sample_df.to_csv(assignment_output, sep='\t')


def dist_nan(u, v):
    valid = ~np.isnan(u) & ~np.isnan(v)
    if valid.sum() == 0:
        return 1
    return cityblock(u[valid], v[valid]) / (2 * valid.sum()) # Manhatten


def euclidean_nan(u, v):
    valid = ~np.isnan(u) & ~np.isnan(v)
    if valid.sum() == 0:
        return 1
    return euclidean(u[valid], v[valid]) / np.sqrt(np.sum(2**2 * valid.sum()))


def main(args):
    df_dna, samples = get_dna_profiles(args.input)
    rel_loci = df_dna.index.drop_duplicates()

    if args.SNParray:
        df_snp, matched_samples = get_snp_array_profiles(args.SNParray, rel_loci)
        label = 'SNP array'
    else:
        df_snp, matched_samples = get_rna_profiles(args.RNA, rel_loci)
        label = 'scRNA-seq'

    df = df_dna.merge(df_snp, left_index=True, right_index=True, how='inner')
    if df.size == 0:
        raise RuntimeError(f'No overlay between snDNA-seq and {label} data!')

    # Check for SNPs with different variant allele
    new_idx = []
    drop_idx = []
    for locus, locus_snps in df_dna.loc[df_snp.index].iterrows():
        locus_df = pd.DataFrame(
            np.reshape(locus_snps.values, (len(matched_samples), 2)),
            columns=['val', 'snp'])

        snps_rel = locus_df[locus_df['val'] > 0.05]['snp'].dropna()
        if snps_rel.nunique() > 1:
            warnings.warn(f'\nSkipping chr{locus[0]}:{locus[1]} due to multiple' \
                f' SNPs: {locus_snps[1::2].dropna().unique()}')
            drop_idx.append(locus)
            continue

        ref, alt = snps_rel.unique()[0].split('>')
        new_idx.append(f'chr{locus[0]}:{locus[1]} {ref}>{alt}')
        for m_sample in matched_samples:
            old_val = df.loc[locus, m_sample]
            try:
                m_ref, m_alt = df.loc[locus, f'snp_{m_sample}'].split('>')
            # Not called/nan
            except AttributeError:
                continue
            # Initial SNP interpretation is fine
            if m_ref == ref and m_alt == alt:
                pass
            # Correct to homozygous snp
            elif m_ref == alt and m_alt == alt:
                if old_val < 0.5:
                    df.at[locus, m_sample] = 1 - old_val
            # Reverse order (dbsnp reference/alternative allele different than reported A/B allele)
            elif m_ref == alt and m_alt == ref:
                df.at[locus, m_sample] = 1 - old_val
            # Wildtype 
            elif m_ref == ref and m_alt == ref:
                if old_val > 0.5:
                    df.at[locus, m_sample] = 1 - old_val
            # No match with dna SNP at all
            elif m_ref != ref and m_ref != alt and m_alt != ref and m_alt != alt:
                df.at[locus, m_sample] = np.nan
            # Different lengths indel
            elif m_ref == ref and len(m_alt) != len(alt):
                df.at[locus, m_sample] = np.nan
            else:
                raise NotImplementedError('Unknown situation:\n' \
                    f'scDNA {ref}->{alt}\nmatched {m_sample}: {m_ref}->{m_alt}')

    df.drop(drop_idx, inplace=True)
    assert len(new_idx) == df.shape[0], \
        'Different SNPs at same locus not implemented yet'

    df = df[sorted(df.columns[::2])]
    df.index = new_idx
    idx_sorted = sorted(df.index,
        key=lambda x: (CHR_ORDER[x.split(':')[0][3:]], x.split(':')[1].split(' ')[0]))
    df = df.loc[idx_sorted]

    names = samples.index.values

    # Remove SNPs that are all wt in the scDNA-seq
    df = df[(df[names] >= 0.1).sum(axis=1) > 0]
    # Remove SNPs that are similar in all DNA clusters
    df = df[(df[names] > 0.95).sum(axis=1) < names.size]
    df = df[(df[names] < 0.05).sum(axis=1) < names.size]
    df = df[((df[names] > 0.45) & (df[names] < 0.55)).sum(axis=1) < names.size]

    dists = []
    for name in names:
        dists.append(np.apply_along_axis(euclidean_nan, 1,
            df[matched_samples.keys()].values.T, df[name].values.T))

    dist_df = pd.DataFrame(dists,index=names, columns=matched_samples)
    dist_df.index.name = 'scDNA-seq'
    dist_df.columns.name = label
    print(f'\nDistance: scDNA-cluster - matched SNPs\n\n{dist_df}\n')
    
    out_base = re.sub(r'(.*)_\d+.filtered_variants.csv$', r'\1', args.input[0])
    if not args.outfile_assignment:
        args.outfile_assignment = f'{out_base}.sample_patient_assigments.tsv'
    assign_clusters(dist_df, samples, matched_samples, args.outfile_assignment)

    if not args.outfile_profiles:
        args.outfile_profiles = f'{out_base}.sample_patient_profiles.png'
    plot_profiles(df.T, args.outfile_profiles)

    if not args.outfile_distance:
        args.outfile_distance = f'{out_base}.sample_patient_distance.png'
    plot_distance(dist_df, args.outfile_distance)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Assigns demultiplexed clusters ' \
        '(<PREFIX>_<CLUSTER>.filtered_variants.csv) to patients based on ' \
        'patient -specific SNPs derived from sc/snRNA-seq or SNP array data.')
    parser.add_argument('-i', '--input', type=str, nargs='+', required=True,
        help='Input scDNAseq variant file(s).')
    parser.add_argument('-r', '--RNA', type=str, nargs='+',
        help='Input scRNAseq SNV profiles.')
    parser.add_argument('-a', '--SNParray', type=str,
        help='Input SNP array file in finalReport.txt (GSGT) format.')
    parser.add_argument('-oa', '--outfile_assignment', type=str,
        help='Output file (format: tsv) for saving the sample assignment.')
    parser.add_argument('-op', '--outfile_profiles', type=str, default='',
        help='Output file (.png|.pdf|...) for profile heatmap.')
    parser.add_argument('-od', '--outfile_distance', type=str, default='',
        help='Output file (.png|.pdf|...) for sample distance heatmap.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.RNA is None and args.SNParray is None:
        raise IOError('No matched SNP data provided.')
    elif args.RNA is not None and args.SNParray is not None:
        warnings.warn('SNP array and scRNA-seq data provided. ' \
            'Using SNP array data.')
    main(args)
