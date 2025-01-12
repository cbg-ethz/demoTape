#!/usr/bin/env python3

import argparse
import os


import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from seaborn import clustermap
import seaborn as sns

EPSILON = np.finfo(np.float64).resolution

FILE_EXT = '.png'
FONTSIZE = 14
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

CHR_ORDER = {str(i): i for i in range(1, 23)}
CHR_ORDER.update({'X':23, 'Y':24})

CL_COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#E4761A', '#2FBF85', '#A4CA55']
DEMOTAPE_COLORS = {
    '0': '#3cb44b', # Green
    '1': '#4363d8', # Blue
    '2': '#ffe119', # Yellow
    '': '#ffffff' # White
} 


def load_cell_annotation(in_file):
    df = pd.read_csv(in_file, index_col=0)
    # DemoTape demultiplexing output
    if df.size == 0:
        df = pd.read_csv(in_file, sep='\t', index_col=0).T
        df['Order'] = df['Order'].astype(int)
        df = df.iloc[df['Order']]
        df.rename({'Cluster': 'assignment'}, axis=1, inplace=True)
        df['assignment'] = df['assignment'].astype(str)
    if 'assignmnet' in df.columns:
        df.rename({'assignmnet': 'assignment'}, axis=1, inplace=True)
    return df


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


def main(args):
    df_in = pd.read_csv(args.input)

    idx_str = df_in['CHR'].astype(str) + ':' + df_in['POS'].astype(str) \
        + ' ' + df_in['REF'] + '>' + df_in['ALT']
    idx_str = idx_str.str.replace('chr', '')

    df = df_in.iloc[:,7:]
    df.index = idx_str

    if not args.output:
        out_base = os.path.splitext(args.input)[0]
    else:
        out_base = os.path.splitext(args.output)[0]

    if args.cell_annotation:
        df_ca = load_cell_annotation(args.cell_annotation)
        assert df.shape[1] == df_ca.shape[0], \
            'Dimension of data and cell annotation inconsistent'
        df = df.loc[:, df_ca.index]
    else:
        df_ca = pd.DataFrame([])
        # Cluster hierarchically
        dist = get_pairwise_dists(df)
        Z = linkage(dist, method='ward')
        cell_order = leaves_list(Z)
        df = df.iloc[:,cell_order]

    plot_heatmap(df, args.output, df_ca)


def plot_heatmap(df, out_file, df_ca):
    try:
        df_plot = df.map(
            lambda x: (int(x.split(':')[1]) + EPSILON) \
                / (sum([int(i) for i in x.split(':')[:2]]) + EPSILON)).T
        mask = df.map(lambda x: x.split(':')[2] == '3').T
    except AttributeError: # Older pandas versions < 2.1.0
        df_plot = df.applymap(
            lambda x: (int(x.split(':')[1]) + EPSILON) \
                / (sum([int(i) for i in x.split(':')[:2]]) + EPSILON)).T
        mask = df.applymap(lambda x: x.split(':')[2] == '3').T

    if df_ca.size > 0:
        # DemoTape assignment
        if df_ca['assignment'].str.contains('+', regex=False).any():
            row_colors = df_ca['assignment'].apply(lambda x: x[0]) \
                .map(DEMOTAPE_COLORS).to_frame(name='Cluster')
            row_colors[''] = df_ca['assignment'] \
                .apply(lambda x: x[-1] if len(x) > 1 else '').map(DEMOTAPE_COLORS)
        # MANTA assignment
        else:
            row_colors = df_ca['assignment'] \
                .map({'tumor': CL_COLORS[0], 'doublets': CL_COLORS[1], \
                    'healthy': CL_COLORS[2]})
    else:
        row_colors = None

    # VAF_cmap = LinearSegmentedColormap.from_list('my_gradient', (
    #     (0.000, (1.000, 1.000, 1.000)),
    #     (0.167, (1.000, 0.882, 0.710)),
    #     (0.333, (1.000, 0.812, 0.525)),
    #     (0.500, (1.000, 0.616, 0.000)),
    #     (0.667, (1.000, 0.765, 0.518)),
    #     (0.833, (1.000, 0.525, 0.494)),
    #     (1.000, (1.000, 0.000, 0.000)))
    # )

    VAF_cmap = LinearSegmentedColormap.from_list('my_gradient', (
        (0.000, (1.000, 1.000, 1.000)),
        (0.500, (1.000, 0.616, 0.000)),
        (1.000, (1.000, 0.000, 0.000)))
    )

    cm = sns.clustermap(
        df_plot,
        mask=mask,
        row_cluster=False,
        row_colors=row_colors,
        col_cluster=False,
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
    if df_plot.shape[0] > 50:
        hm.set_yticks([])

    hm.set_xlabel('SNPs')
    hm.set_xticklabels(hm.get_xticklabels(), fontsize=FONTSIZE/2, va='top')

    cm.ax_cbar.set_title('VAF', fontsize=FONTSIZE)

    if out_file:
        cm.fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        help='Variant file from mosaic preprocessing.'),
    parser.add_argument('-ca', '--cell_annotation', type=str, default='',
        help='Cell order (barcodes) in csv format.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file. Default = stdout.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)