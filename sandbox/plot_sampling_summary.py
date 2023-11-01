#!/usr/bin/env python3

import argparse
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import seaborn as sns


COLORS = {
    'BnpC': '#8534B8',
    'SCG (8% dbt)': '#90E83D',
    'SCG (30% dbt)': '#4D9D00',
    'COMPASS (8% dbt)': '#FF9543',
    'COMPASS (30% dbt)': '#D75D00',
    '': '#FFFFFF'
}

SID2PAPER = {'SID3841': 'S1', 'SID3867': 'S2', 'SID5962': 'S3'}
MAX_CELLS = {'S1': 3732, 'S2': 7911, 'S3': 3677}
PANELS = [f'{i}' for i in map(chr, range(65, 91))] # 97, 123 for small


FIG_SIZE = (12, 5)
FILE_EXT = 'png'
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



def plot_cluster_no(df, ax):
    n_alg = df['algorithm'].nunique()

    bp = sns.boxplot(
        data=df,
        x='cells',
        y='#Cluster',
        hue='algorithm',
        ax=ax,
        fliersize=2,
        showfliers=False,
        linewidth=1,
        palette=COLORS
    )
    sns.stripplot(
        data=df,
        x='cells',
        y='#Cluster',
        hue='algorithm',
        ax=ax,
        linewidth=1,
        jitter=0.15,
        alpha=.75,
        size=4,
        dodge=True,
        palette=COLORS
    )
    ax.set_ylabel('#Clusters')
    ax.set_ylim((0, df['#Cluster'].max() + 1))
    ax.set_yticks(range(0, df['#Cluster'].max() + 1, 2))

    ax.set_xlabel('Sample size [cells]')
    ax.set_xticks(range(df['cells'].nunique()))
    x_tick_labels = [str(j) if i % 2 != 0 else '' \
        for i, j in enumerate(np.unique(df['cells']))]
    ax.set_xticklabels(x_tick_labels)

    labels = ['', 'SCG (8% dbt)', 'SCG (30% dbt)', 'BnpC',
        'COMPASS (8% dbt)', 'COMPASS (30% dbt)', ]
    ax.legend(
        [mpatches.Patch(color=COLORS[i]) for i in labels],
        labels,
        title=r'$\bf{Algorithm}$',
        bbox_to_anchor=(0., 1.02, 1., .102),
        loc='lower left',
        ncol=2,
        mode="expand",
        borderaxespad=0.
    )


def plot_geno_dist(df, ax):
    n_alg = df['algorithm'].nunique()
    bp = sns.boxplot(
        data=df,
        x='cells',
        y='Genotyping Sim.',
        hue='algorithm',
        ax=ax,
        fliersize=2,
        showfliers=False,
        linewidth=1,
        palette=COLORS
    )
    sns.stripplot(
        data=df,
        x='cells',
        y='Genotyping Sim.',
        hue='algorithm',
        ax=ax,
        linewidth=1,
        jitter=0.15,
        alpha=.75,
        size=6,
        dodge=True,
        palette=COLORS
    )

    ax.set_ylim((0.75, 1))
    ax.set_yticklabels([0, 0.80, 0.85, 0.90, 0.95, 1.00])

    x_tick_labels = [str(j) if i % 2 != 0 else '' \
        for i, j in enumerate(np.unique(df['cells']))]
    ax.set_xticklabels(x_tick_labels)

    # Plot linebreak
    d = 0.02
    y_br = 0.06
    br_width = 0.02
    ax.plot((-d, +d), (y_br-d, y_br+d),
        c='k', clip_on=False, transform=ax.transAxes, lw=2)
    ax.plot((-d, +d), (y_br-d+br_width, y_br+d+br_width),
        c='k', clip_on=False, transform=ax.transAxes, lw=2)

    ax.set_ylabel(r'Genotyping consistency $\left( 1 - \frac{| G - G^{\prime}|}{norm.} \right)$')
    ax.set_xlabel('Sample size [cells]')

    ax.get_legend().remove()


def plot_coclustering_dist(df, ax):
    n_alg = df['algorithm'].nunique()
    bp = sns.boxplot(
        data=df,
        x='cells',
        y='Clustering Sim.',
        hue='algorithm',
        ax=ax,
        fliersize=2,
        showfliers=False,
        linewidth=1,
        palette=COLORS
    )
    sns.stripplot(
        data=df,
        x='cells',
        y='Clustering Sim.',
        hue='algorithm',
        ax=ax,
        linewidth=1,
        jitter=0.15,
        alpha=.75,
        size=6,
        dodge=True,
        palette=COLORS
    )

    ax.set_ylim((0.2, 1))
    ax.set_yticklabels([0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    x_tick_labels = [str(j) if i % 2 != 0 else '' \
        for i, j in enumerate(np.unique(df['cells']))]
    ax.set_xticklabels(x_tick_labels)

    # Plot linebreak
    d = 0.02
    y_br = 0.06
    br_width = 0.02
    ax.plot((-d, +d), (y_br-d, y_br+d),
        c='k', clip_on=False, transform=ax.transAxes, lw=2)
    ax.plot((-d, +d), (y_br-d+br_width, y_br+d+br_width),
        c='k', clip_on=False, transform=ax.transAxes, lw=2)

    ax.set_ylabel(r'Clustering consistency $\left( 1 - \frac{| M - M^{\prime}|}{norm.} \right)$')
    ax.set_xlabel('Sample size [cells]')

    ax.get_legend().remove()


def plot_metric(metric, df, ax):
    n_alg = df['algorithm'].nunique()
    bp = sns.boxplot(
        data=df,
        x='cells',
        y=metric,
        hue='algorithm',
        ax=ax,
        fliersize=2,
        showfliers=False,
        linewidth=1,
        palette=COLORS
    )
    sns.stripplot(
        data=df,
        x='cells',
        y=metric,
        hue='algorithm',
        ax=ax,
        linewidth=1,
        jitter=0.15,
        alpha=.75,
        size=6,
        dodge=True,
        palette=COLORS
    )

    x_tick_labels = [str(j) if i % 2 != 0 else '' \
        for i, j in enumerate(np.unique(df['cells']))]
    ax.set_xticklabels(x_tick_labels)

    # ax.set_ylim((0.0, 1))
    # ax.set_yticklabels([0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    # Plot linebreak
    # d = 0.02
    # y_br = 0.06
    # br_width = 0.02
    # ax.plot((-d, +d), (y_br-d, y_br+d),
    #     c='k', clip_on=False, transform=ax.transAxes, lw=2)
    # ax.plot((-d, +d), (y_br-d+br_width, y_br+d+br_width),
    #     c='k', clip_on=False, transform=ax.transAxes, lw=2)

    if metric == 'ARI':
        ax.set_ylabel(r'Adjusted Rand Index')
    elif metric == 'Mutual info':
        ax.set_ylabel(r'Adjusted Mutual Information')
    ax.set_xlabel('Sample size [cells]')

    ax.get_legend().remove()



def main(args):
    col_no = len(args.input)
    first_df =  pd.read_csv(args.input[0], sep='\t')
    min_cells_all = first_df['min. #cells per cluster'].unique()
    metrics = ['#Cluster', 'ARI']
    row_no = len(metrics)


    for min_cells in min_cells_all:
        fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
            figsize=(FIG_SIZE[0]*col_no, FIG_SIZE[0]*row_no))
        axes = np.reshape(axes, (row_no, col_no))

        for col, in_file in enumerate(args.input):
            df = pd.read_csv(in_file, sep='\t')
            df.drop_duplicates(inplace=True)
            try:
                df = df[df['min. #cells per cluster'] == min_cells]
            except:
                pass

            if 'doublet rate' in df:
                df['algorithm'] = df.apply(lambda x:
                    f'{x["algorithm"]} ({x["doublet rate"] * 100:.0f}% dbt)' \
                        if x["doublet rate"] > 0 else f'{x["algorithm"]}' ,axis=1)

            max_c = df['cells'].unique().max()
            name = os.path.basename(in_file).split('_')[0]
            if name in SID2PAPER:
                name = SID2PAPER[name]
                df['cells'].replace(max_c, MAX_CELLS[name], inplace=True)
                max_c = MAX_CELLS[name]

            plot_cluster_no(df, axes[0][col])
            add_panel_id(axes[0][col], PANELS[0*col_no + col])
            axes[0][col].set_title(f'{name}\n\n\n\n', fontsize='xx-large', fontweight='bold')

            plot_metric('ARI', df, axes[1][col])
            add_panel_id(axes[1][col], PANELS[1*col_no + col])

            # plot_coclustering_dist(df, axes[1][col])
            # add_panel_id(axes[1][col], PANELS[1*col_no + col])
            # plot_geno_dist(df, axes[2][col])
            # add_panel_id(axes[2][col], PANELS[2*col_no + col])
            # plot_metric('Mutual info', df, axes[-1][col])
            # add_panel_id(axes[-1][col], PANELS[2*col_no + col])


        fig.tight_layout()
        if args.output:
            out_file = f'{{}}.min{min_cells}perCl{{}}' \
                .format(*os.path.splitext(args.output))
        else:
            out_file = os.path.join(os.path.dirname(args.input[0]),
                f'sampling.min{min_cells}perCl.{FILE_EXT}')

        print(f'Writing file to: {out_file}')
        fig.savefig(out_file, dpi=DPI)


def add_panel_id(ax, label):
    ax.text(-0.1, 1.05, label, transform=ax.transAxes, fontsize='xx-large',
        va='bottom', fontweight='bold', ha='center')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Sampling summary tsv file.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for result png. default = <DIR:INPUT[0]>/sampling.png')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)

