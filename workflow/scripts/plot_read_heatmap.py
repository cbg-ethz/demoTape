#!/usr/bin/env python3

import argparse
from itertools import cycle
import os
import re

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import f_oneway
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
        'axes.titlesize': 'medium',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 10
})
plt.rcParams['xtick.major.size'] = 1
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['xtick.bottom'] = False

CHR_ORDER = {str(i): i for i in range(1, 23)}
CHR_ORDER.update({'X':23, 'Y':24})

COLOR_GENE_BUF = cycle(['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
    '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928'])
COLOR_CHROM_BUF = cycle(['#f4f4f4','#c3c4c3'])
CL_COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#E4761A' , '#2FBF85', '#A4CA55']


def calc_gini(x):
    total = 0
    for i, xi in enumerate(x[:-1], 1):
        total += np.sum(np.abs(xi - x[i:]))
    return (total + EPSILON) / (len(x)**2 * np.mean(x) + EPSILON)


def load_cell_annotation(in_file):
    df = pd.read_csv(in_file, index_col=0)
    if 'assignmnet' in df.columns:
        df.rename({'assignmnet': 'assignment'}, axis=1, inplace=True)
    return df


def load_panel(in_file):
    df = pd.read_csv(in_file, comment='#', sep='\t', header=None, index_col=-1)
    df[0] = df[0].str.replace('chr', '')
    df['locus'] = df[0].astype(str) + ':' + df[1].astype(str) \
            + '-' + df[2].astype(str)
    return df


def get_pairwise_dists(df):
    dp = ((df * 1e6).round().astype(int) + 1).values
    dp_fct = dp * np.log(dp)
    norm_const = np.arange(0, dp.max().max() * 2 + 1) * np.log(2)

    dist = []
    for i in np.arange(dp.shape[0] - 1):
        valid = (dp[i] > 1) & (dp[i+1:] > 1)
        dp_tot = dp[i] + dp[i+1:]
        l12 = dp_tot / 2
        logl = dp_fct[i] + dp_fct[i+1:] \
            - (dp_tot) * np.log(l12) \
            + 2 * l12 - dp_tot
        norm = norm_const[dp_tot]
        dist.append(np.sum(np.where(valid, logl / norm, 0), axis=1) \
            / valid.sum(axis=1))
    return np.concatenate(dist)



def main(args):
    df = pd.read_csv(args.input, sep='\t', header=0, index_col=0)
    if args.sample_annotation:
        df_sa = pd.read_csv(args.sample_annotation, sep='\t', header=0, index_col=0)
        try:
            sa_id = int(re.search(r'_(\d)\.', args.output).group(1))
        except AttributeError: # Original, still multiplexed, sample
            sa_id = -1
        else:
            sa_cells = df_sa.loc[:, df_sa.loc['Cluster'] == sa_id].columns.values
            df = df.loc[sa_cells]
    df.sort_index(inplace=True)

    if args.panel:
        panel = load_panel(args.panel)
        # Order amplicons based on their chromosomal positions
        ampl_order = df.rename(panel['locus'].to_dict(), axis=1).columns.values
        ampl_idx_sorted = sorted(range(df.shape[1]),
            key=lambda i:(CHR_ORDER[ampl_order[i].split(':')[0]], i))
        df = df[df.columns[ampl_idx_sorted]]
    else:
        panel = pd.DataFrame([])

    lib_depth = np.log10(df.sum(axis=1)).to_frame(name='Library\nsize [log]')
    # df.clip(upper=1000, axis=1, inplace=True)
    # High expression amplicons
    high_exp = df.mean() > 2 * df.mean().mean()

    # Noisy amplicons (tapestri)
    gini = df.apply(lambda i: calc_gini(i.values))
    gini_threshold = gini.mean() + (2 * gini.std())
    noisy = gini > gini_threshold

    # low performing amplicons (tapestri)
    threshold = 0.2 * df.mean().mean()
    low_amplicons = df.mean() < threshold
    
    # Remove amplicons that are covered in less than 20% of all cells
    df = df.loc[:, (df == 0).sum() / df.shape[0] <= 0.8]
    # Remove amplicons with <= 1 read / cell
    low_cov = (df.sum() / df.shape[0]) < 1

    # normalize read counts per cell
    df = df.apply(lambda x: (x / x.sum()), axis=1)
    # drop bad amplicons
    df = df.loc[:,~(noisy | low_amplicons | low_cov)]
    panel = panel.loc[~(noisy | low_amplicons | low_cov)]

    # Remove outliers: clip data to 10% and 90% quantile
    df.clip(lower=None, upper=df.quantile(0.9), axis=1, inplace=True)

    if args.cell_annotation:
        df_ca = load_cell_annotation(args.cell_annotation)
        if args.sample_annotation and sa_id >= 0:
            df_ca = df_ca.loc[sa_cells]
        df = df.loc[df_ca.index]
        lib_depth = lib_depth.loc[df_ca.index]
        lib_depth['assignment'] = df_ca['assignment'].values
    else:
        # Cluster hierarchically
        dist = get_pairwise_dists(df)
        Z = linkage(dist, method='ward')
        cell_order = leaves_list(Z)
        df = df.iloc[cell_order]

    # Normalize per amplicon
    df = df / df.mean()

    # Normalize such that avg. healthy cell depth = 2
    if args.cell_annotation:
        tumor_cl = df_ca[df_ca['assignment'] == 'tumor'].iloc[0]['cluster']
        healthy_cl = df_ca[df_ca['assignment'] == 'healthy'].iloc[0]['cluster']
        
        tumor_cells = df_ca['cluster'] == tumor_cl
        healthy_cells = df_ca['cluster'] == healthy_cl

        df = df.apply(lambda x: x / x[healthy_cells].mean() * 2, axis=0)
    else:
        tumor_cells = df.index

    # Output CNV profiles for tools like cloneAlign or ScaTrex
    CNV_profiles = df.loc[tumor_cells].T
    CNV_profiles['gene'] = panel.loc[CNV_profiles.index][3]
    
    CNV_profiles = CNV_profiles.groupby('gene').mean().mean(axis=1) \
        .round().astype(int).to_frame(name='tumor')
    CNV_profiles['healthy'] = 2
    CNV_profiles.drop('.', inplace=True, errors='ignore') # Drop amplicons not covering a gene
    
    prefix = os.path.basename(args.input).split('.')[0]
    if args.sample_annotation and sa_id >= 0:
        prefix += f'_{sa_id}'
    if args.output:
        out_dir = os.path.dirname(args.output)
    else:
        out_dir = os.path.dirname(args.input)
    profiles_out_file = os.path.join(out_dir, f'{prefix}.reads_CNVclones.csv')
    print(f'Writing CNV profiles to: {profiles_out_file}')
    CNV_profiles.to_csv(profiles_out_file)

    # Normalize such that avg. cell depth = 2
    if not args.cell_annotation:
        df = df.apply(lambda x: x / x.mean() * 2, axis=0)

    plot_heatmap(df.round(1), args.output, panel, lib_depth, args.show_ampl)


def plot_heatmap(df, out_file, panel, lib_depth, show_ampl=False):
    if panel.size > 0:
        chr_colors = []
        c_color_map = {}
        gene_colors = []
        g_color_map = {}

        for ampl in df.columns:
            gene = panel.loc[ampl][3]
            chrom = panel.loc[ampl][0]
            if gene not in g_color_map:
                g_color_map[gene] = next(COLOR_GENE_BUF)
            if chrom not in c_color_map:
                c_color_map[chrom] = next(COLOR_CHROM_BUF)
            gene_colors.append(g_color_map[gene])
            chr_colors.append(c_color_map[chrom])

        if show_ampl:
            col_data =  [chr_colors, chr_colors]
            idx = ['Chr', '']
        else:
            col_data = [gene_colors, chr_colors]
            idx = ['Gene', 'Chr']
            
        col_colors = pd.DataFrame(col_data, columns=df.columns, index=idx).T
    else:
        col_colors = None

    norm = Normalize(vmin=lib_depth['Library\nsize [log]'].min(),
        vmax=lib_depth['Library\nsize [log]'].max(), clip=True)
    mapper = ScalarMappable(norm=norm, cmap='afmhot')
    row_colors = lib_depth['Library\nsize [log]'].apply(mapper.to_rgba).to_frame()
    if 'assignment' in lib_depth:
        row_colors['assignment'] = lib_depth['assignment'] \
            .map({'tumor': CL_COLORS[0], 'doublets': CL_COLORS[1], \
                'healthy': CL_COLORS[2]})

    CNV_cmap = LinearSegmentedColormap.from_list('my_gradient', (
        (0.000, (0.188, 0.400, 0.776)), # Dark Blue
        (0.250, (0.106, 0.376, 1.000)), # Light Blue
        (0.450, (1.000, 1.000, 1.000)), # White
        (0.563, (1.000, 1.000, 1.000)), # White
        (0.750, (1.000, 0.000, 0.000)), # Red
        (1.000, (0.729, 0.000, 0.000))) # Dark Red
    )

    cm = clustermap(
        df,
        row_cluster=False,
        row_colors=row_colors,
        col_cluster=False,
        col_colors=col_colors,
        cmap=CNV_cmap,
        vmin=0, center=2, vmax=6,
        figsize=(25, 10),
        cbar_kws={'shrink': 0.5, 'ticks': [0, 1, 2, 3, 4, 5, 6]},
        cbar_pos=(0.15, 0.8, 0.01, 0.075),
        xticklabels=True
    )

    hm = cm.ax_heatmap

    hm.set_ylabel('Cells')
    if df.shape[0] > 50:
        hm.set_yticks([])

    if show_ampl:
        x_labels = [panel.loc[i.get_text(), 'locus'] for i in hm.get_xticklabels()]
        x_labels = [i.get_text() for i in hm.get_xticklabels()]
    else:
        x_labels = ['{} (chr{})'.format(*panel.loc[i.get_text(), [3, 0]]) \
            for i in hm.get_xticklabels()]
    hm.set_xlabel('Amplicons')
    hm.set_xticklabels(x_labels, fontsize=2, va='top')

    # Remove position of gene and chrom column colors
    box_cc = cm.ax_col_colors.get_position()
    cm.ax_col_colors.set_position([box_cc.x0,
        hm.get_position().min[1] - 0.7 * box_cc.height - 0.005,
        box_cc.width, 0.7 * box_cc.height])

    cm.ax_cbar.set_title('Norm. read depth', fontsize=FONTSIZE/2)
    for i in cm.ax_cbar.get_yticklabels():
        i.set_fontsize('x-small')

    if out_file:
        print(f'Writing read heatmap to: {out_file}')
        cm.savefig(out_file, dpi=DPI)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
        help='Reads per barcode distribution from Tapestri processing (.tsv).'),
    parser.add_argument('-p', '--panel', type=str, 
        default='resources/Tapestri-Designer-results-4387/4387_annotated.bed',
        help='Tapestri panel bed file')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for summary tsv. default = <DIR:INPUT>/read_ampl_depth.png')
    parser.add_argument('-sa', '--sample_annotation', type=str, default='',
        help='Cell-to-sample annotation for demultiplexed samples.')
    parser.add_argument('-ca', '--cell_annotation', type=str, default='',
        help='Cell order (barcodes) in csv format.')
    parser.add_argument('-a', '--show_ampl', action='store_true',
        help='Show loci (chr:pos) instead of gene names.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)