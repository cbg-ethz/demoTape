#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from seaborn import clustermap
import seaborn as sns

EPSILON = np.finfo(np.float64).resolution

FILE_EXT = '.png'
FONTSIZE = 16
DPI = 300
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'lines.linewidth': 2,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'large',
        'xtick.labelsize': 'medium',
        'xtick.major.size': 1,
        'xtick.major.width': 0.5,
        'xtick.bottom': False,
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 10
})

CHR_ORDER = dict({str(i): i for i in range(1, 23, 1)}, **{'X': 23, 'Y': 24})

PATIENT_COLORS = {
    'S3': '#5FCA6C',
    'MS1:S3': '#0d821b',
    'S1': '#6580E2',
    'MS1:S1': '#0e2fa7',
    'S2': '#ffe119',
    'MS1:S2': '#a69100',
    'MS1': '#6F41D9',
}
CALLED_COLORS = {0: '#FFF0F0', 1: '#81A281'}


def get_pairwise_dists(df):
    ref = df.T.map(lambda x: int(x.split(':')[0])).values
    alt = df.T.map(lambda x: int(x.split(':')[1])).values
    dp = ref + alt
    VAF = np.clip((alt + EPSILON) / (dp + EPSILON), EPSILON, 1 - EPSILON)
    RAF = 1 - VAF

    norm_const = np.insert(
        np.arange(1, dp.max().max() * 2 + 1) \
            * np.log(np.arange(1, dp.max().max() * 2 + 1)),
        0, np.nan)

    dist = []
    for i in np.arange(VAF.shape[0] - 1):
        valid = (dp[i] > 0) & (dp[i+1:] > 0)
        dp_total = dp[i] + dp[i+1:]
        p12 = np.clip((alt[i] + alt[i+1:] + EPSILON) / (dp_total + EPSILON),
            EPSILON, 1 - EPSILON)
        p12_inv = 1 - p12
        logl = alt[i] * np.log(VAF[i] / p12) \
            + ref[i] * np.log(RAF[i] / p12_inv) \
            + alt[i+1:] * np.log(VAF[i+1:] / p12) \
            + ref[i+1:] * np.log(RAF[i+1:] / p12_inv)

        norm = norm_const[dp_total] \
            - norm_const[dp[i]] \
            - norm_const[dp[i+1:]]

        dist.append(np.nansum(logl / norm, axis=1) / valid.sum(axis=1))
    return np.concatenate(dist)


def get_idx(df, short=False):
    idx = df['CHR'].astype(str) + ':' + df['POS'].astype(str) 
    if short:
        return idx.str.replace('chr', '')
    return idx.str.replace('chr', '') + ' ' + df['REF'] + '>' + df['ALT']


def main(args):
    s_ca = pd.Series([], name='Sample')
    for in_file in args.input_whitelist:
        df_in = pd.read_csv(in_file)
        idx_str = get_idx(df_in)
        df_in.index = idx_str

        # Order cells per sample based on dendrogram
        df_in = df_in.iloc[:,7:]
        dist = get_pairwise_dists(df_in)
        Z = linkage(dist)
        cell_order = leaves_list(Z)
        df_in = df_in[df_in.columns[cell_order]]
        
        sample = os.path.basename(in_file).split('.')[1]
        if 'demulti' in os.path.basename(in_file):
            sample = f'MS1:{sample}'
        df_in.columns += f'_{sample}'

        new_s_ca = pd.Series(sample, index=df_in.columns, name='Sample')
        s_ca = pd.concat([s_ca, new_s_ca])

        try:
            df = df.merge(df_in, left_index=True, right_index=True, how='outer')
            # Check for duplicated loci with different alternative alleles and merge
            df['loci'] = [i.split(' ')[0] for i in df.index]
            dbl_loci = np.argwhere(df['loci'].duplicated(keep=False) == True).ravel()
            for dbl_locus1 in dbl_loci[::2]:
                dbl_locus2 = dbl_locus1 + 1
                locus_comb = df.iloc[dbl_locus1].fillna('') \
                    + df.iloc[dbl_locus2].fillna('')
                df.iloc[dbl_locus1] = locus_comb
                df.iloc[dbl_locus2] = locus_comb
            df.drop_duplicates(inplace=True)
            df.drop('loci', axis=1, inplace=True)
        except NameError:
            df = df_in

    idx_sorted = sorted(df.index, key=lambda x: (CHR_ORDER[x.split(':')[0]], x))
    df = df.loc[idx_sorted]

    # Get df indicating which SNPs were called in the demultiplexed and single sample
    called_snps = pd.DataFrame(0,
        index=[i.split(' ')[0] for i in df.index], columns=s_ca.unique())
    for orig_file in args.input_original:
        df_orig =  pd.read_csv(orig_file)
        idx_str = get_idx(df_orig, short=True)
        # demultiplexed file
        if re.search(r'_\d\.(?=relevant|filtered)', orig_file):
            sample = [i for i in called_snps.columns if i.startswith('MS1')][0]
        # single file
        else:
            sample = [i for i in called_snps.columns if not i.startswith('MS1')][0]
        called_snps.loc[idx_str, sample] = 1

    called_snps.index = df.index
    called_snps.columns = 'called in ' + called_snps.columns

    plot_heatmap(df, s_ca, called_snps, args.output)


def plot_heatmap(df, s_ca, df_called, out_file):
    ref = df.T.map(lambda x: int(x.split(':')[0])).values
    alt = df.T.map(lambda x: int(x.split(':')[1])).values

    dp = ref + alt
    VAF = np.clip((alt + EPSILON) / (dp + EPSILON), EPSILON, 1 - EPSILON)
    df_plot = pd.DataFrame(VAF, index=s_ca.index, columns=df_called.index)

    mask = dp == 0
    col_colors = df_called.replace(CALLED_COLORS)
    row_colors = s_ca.replace(PATIENT_COLORS)

    VAF_cmap = LinearSegmentedColormap.from_list('my_gradient', (
        (0.000, (1.000, 1.000, 1.000)),
        (0.500, (1.000, 0.616, 0.000)),
        (1.000, (1.000, 0.000, 0.000)))
    )

    cm = sns.clustermap(
        df_plot,
        mask=mask,
        row_cluster=False,
        row_colors=row_colors.to_frame(),
        col_cluster=False,
        col_colors=col_colors,
        cmap=VAF_cmap,
        vmin=0, center=0.5, vmax=1,
        figsize=(25, 10),
        cbar_kws={'shrink': 0.5, 'ticks': [0, 0.5, 1]},
        cbar_pos=(0.15, 0.8, 0.01, 0.075),
        xticklabels=True
    )

    hm = cm.ax_heatmap
    hm.set_facecolor('#636161')

    hm.set_ylabel('Cells')
    if VAF.shape[0] > 50:
        hm.set_yticks([])

    hm.set_xlabel('SNPs')
    hm.set_xticklabels(hm.get_xticklabels(), fontsize=int(FONTSIZE * 0.75), va='top')

    cm.ax_cbar.set_title('VAF', fontsize=FONTSIZE)
    # plt.show()
    # exit()

    if out_file:
        cm.fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-iwl', '--input_whitelist', type=str, nargs='+',
        help='Variant files from mosaic preprocessing after whitelisting.'),
    parser.add_argument('-io', '--input_original', type=str, nargs='+',
        help='Variant files from mosaic preprocessing before whitelisting.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file. Default = stdout.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)