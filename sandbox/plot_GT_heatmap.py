#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from seaborn import clustermap
import seaborn as sns

CHR_ORDER = {str(i): i for i in range(1, 23)}
CHR_ORDER.update({'X':23, 'Y':24})

COLORS = {0: '#e41a1c', 1: '#377eb8', 2: '#4daf4a', 3: '#E4761A'}
COLORS0 = {1: '#e41a1c', 0: '#FF7575'} # red
COLORS1 = {1: '#377eb8', 0: '#84B5DE'} # blue
COLORS2 = {1: '#4daf4a', 0: '#AAE3A7'} # green
COLORS = COLORS0


FILE_EXT = '.png'
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
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.bottom'] = True


def get_whitelist_loci(wl_file):
    loci = []
    df = pd.read_csv(wl_file, sep='\t')

    if df.size == 0:
        if df.shape[1] > 1:
            for snv in df.columns:
                loci.append('_'.join(snv.split('_')[:2]))
        else:
            for snv in df.columns.values[0].split(';'):
                loci.append('_'.join(snv.split('_')[:2]))
    elif df.shape[1] == 1:
        df = pd.read_csv(wl_file, sep=',')
        ids = df['CHR'].astype(str) + '_' + df['POS'].astype(str)
        loci.extend(ids)
    else:
        for _, (sample, snvs) in df.iterrows():
            if '+' in sample:
                continue
            for snv in snvs.split(';'):
                loci.append('_'.join(snv.split('_')[:2]))
    return loci



def get_row_colors(gt, cluster_assignment):
    if cluster_assignment:
        cluster_rows = False
        assignment = pd.read_csv(cluster_assignment,
            index_col=0, dtype=str, sep='\t')
        assign_vec = assignment.values.flatten()
        cell_order = []
        row_colors = [[], []]
        #for i, cl in enumerate(np.unique(assign_vec)):
        for i, cl in enumerate(['1', '1+2', '0+1', '2', '0+2', '0']):
            cell_ids = np.argwhere(assign_vec == cl).flatten()
            cell_order.extend(cell_ids)
            if '+' in cl:
                s1, s2 = [int(j) for j in cl.split('+')]
                if s1 > s2:
                    row_colors[0].extend(np.full(cell_ids.size, COLORS[s2]))
                    row_colors[1].extend(np.full(cell_ids.size, COLORS[s1]))
                else:
                    row_colors[0].extend(np.full(cell_ids.size, COLORS[s1]))
                    row_colors[1].extend(np.full(cell_ids.size, COLORS[s2]))
            else:
                row_colors[0].extend(np.full(cell_ids.size, COLORS[int(cl)]))
                row_colors[1].extend(np.full(cell_ids.size, '#FFFFFF'))
        gt = gt.iloc[:, cell_order]
    else:
        cluster_rows = False
        file_ids = [int(i[-1]) for i in gt.columns]
        cl_no = np.unique(file_ids).size
        # Get row color colors based on assignment
        if cl_no > 1:
            row_colors = [*map(COLORS.get, file_ids)]
        else:
            row_colors = None
    return cluster_rows, row_colors, gt


def main(args):
    for i, in_file in enumerate(args.input):
        df = pd.read_csv(in_file, sep=',')
        if 'CHR' in df.columns:
            too_long = df['REF'].str.len() > 10
            df['REF'] = df['REF'].str[:10]
            df.loc[too_long, 'REF'] = df[too_long]['REF'] + '...'

            ids = df.apply(lambda x: f'{x["CHR"]}_{x["POS"]}_{x["REF"]}_{x["ALT"]}',
                axis=1)
            df.set_index(ids, inplace=True)
            df.columns += f'.{i}'
            gt_new = df.iloc[:,7:].applymap(lambda x: float(x.split(':')[-1]))
        else:
            df = pd.read_csv(in_file, sep=',',index_col=0).T
            df.columns = [f'{j}.{i}' for j in df.columns]
            gt_new = df

        try:
            gt = gt.merge(gt_new, how='outer', left_index=True, right_index=True)
        except NameError:
            gt = gt_new

    idx_sorted = sorted(gt.index,
        key=lambda i:(CHR_ORDER[i.split('_')[0]], int(i.split('_')[1])))
    gt = gt.loc[idx_sorted]

    # unique_sites = gt[gt.isna().sum(axis=1) > 0]
    # unique_sites = unique_sites.replace(3, np.nan)
    # cl_unique = unique_sites.mean(axis=1) > 0.1
    # import pdb; pdb.set_trace()

    gt.fillna(-1, inplace=True)
    gt.replace(3, -1, inplace=True)

    if args.whitelist:
        wl_loci = get_whitelist_loci(args.whitelist)
        gt['locus'] = ['_'.join(i.split('_')[:2]) for i in gt.index]
        gt = gt[gt['locus'].isin(wl_loci)]
        gt = gt.drop('locus', axis=1)

    if args.filter_wildtype:
        wt_frac = (gt == 0).sum(axis=1) / (gt >= 0).sum(axis=1)
        gt = gt[wt_frac < 0.95]

    # Remove SNVs with >95% missing
    # gt = gt[(gt == -1).mean(axis=1) < 0.95]

    cluster_rows, row_colors, gt = get_row_colors(gt, args.cluster_assignment)

    plot_heatmap(gt.T, args.output, cluster_rows, row_colors)


def plot_heatmap(gt, out_file, cluster_rows, row_colors):
    myColors = ('#EAEAEA', '#fed976', '#fc4e2a', '#800026')
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    cm = clustermap(
        gt,
        row_cluster=cluster_rows,
        col_colors=None,
        col_cluster=True,
        row_colors=row_colors,
        vmin=-1, vmax=2,
        cmap=cmap,
        dendrogram_ratio=(0.05, 0.3),
        figsize=(25, 10),
        cbar_kws={'shrink': 0.5, 'drawedges': True},
        cbar_pos=(0, 0, 0.0, 0.0),
        xticklabels=True,
        tree_kws=dict(linewidths=2)
    )
    cm.ax_heatmap.set_facecolor('#EAEAEA')

    cm.ax_heatmap.set_ylabel('Cells')
    if gt.shape[0] > 50:
        cm.ax_heatmap.set_yticks([])

    if gt.shape[1] > 25:
        cm.ax_col_dendrogram.set_visible(False) # Turn col heatmap off
        label_font = 15
    else:
        label_font = FONTSIZE

    cm.ax_row_dendrogram.set_visible(False)

    cm.ax_heatmap.set_xlabel('SNVs')
    labels_pretty = ['chr{}:{} {}>{}'.format(*i.get_text().split('_')) \
        for i in cm.ax_heatmap.get_xticklabels()]

    cm.ax_heatmap.set_xticklabels(labels_pretty,
        rotation=45, fontsize=label_font, ha='right', va='top')

    cm.ax_cbar.set_title('Genotype')
    colorbar = cm.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([-0.65, 0.15, 0.85, 1.65])
    colorbar.set_ticklabels([r' $-$', '0|0', '0|1', '1|1'])


    cm.fig.tight_layout()
    if out_file:
        # if not args.output_plot.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
        #     args.output_plot += '.png'
        cm.fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Variant files from mosaic preprocessing.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for summary tsv. default = <DIR:INPUT[0]>/summary.tsv')
    parser.add_argument('-wl', '--whitelist', type=str, default='',
        help='Whitelist containint SNVs for plotting. Default = None.')
    parser.add_argument('-cl', '--cluster_assignment', type=str, default='',
        help='Cluster assignment for multiplex variant file.')
    parser.add_argument('-fwt', '--filter_wildtype', action='store_true',
        help='Remove SNPs with >95 wildtype genotype.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)