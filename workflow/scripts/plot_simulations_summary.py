#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context('talk')
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.transforms as mtransforms

MR_MAP = {
    '0.33-0.33-0.33': 0,
    '0.2-0.4-0.4': 1, '0.4-0.2-0.4': 1, '0.4-0.4-0.2': 1,
    '0.2-0.3-0.5': 2, '0.3-0.2-0.5': 2, '0.5-0.3-0.2': 2, '0.5-0.2-0.3': 2,
    '0.2-0.5-0.3': 2, '0.3-0.5-0.2': 2,
    '0.1-0.3-0.6': 3, '0.3-0.1-0.6': 3, '0.6-0.3-0.1': 3, '0.6-0.1-0.3': 3,
    '0.1-0.6-0.3':3, '0.3-0.6-0.1': 3
}
MR_MAP_reverse = {0: '33:33:33%', 1: '20:40:40%', 2: '20:30:50%', 3: '10:30:60%'}
ALG_MAP = {
    'scSplit': 'scSplit',
    'Euclidean': 'demoTape',
    'souporcell': 'souporcell',
    'demoTape': 'demoTape',
    'vireo': 'Vireo'
}
PANELS = [f'{i}' for i in map(chr, range(65, 91))] # 97, 123 for small
COLORS = {
    'demoTape': '#F24C3D',
    'scSplit': '#F2BE22',
    'souporcell': '#22A699',
    'vireo': '#325FB0',
}
METRIC = 'v_measure'

MARGINS = {
    'left': 0.05,
    'bottom': 0.1,
    'right': 1,
    'top': 0.9,
    'wspace': 0.3,
    'hspace': 0.3 
}


FONTSIZE = 50
DPI = 300
sns.set_style('whitegrid')
sns.set_context('paper',
    rc={'xtick.major.size': 2,
        'ytick.major.size': 2,
        'lines.linewidth': 2,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'large',
        'axes.titlesize': 'large',
        'xtick.labelsize': 'large',
        'ytick.labelsize': 'large',
        'legend.fontsize': 'large',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 50,
        'lines.linewidth': 1,
        'xtick.major.size':  6,
        'ytick.major.size':  6,
})


def main(in_file, out_file):
    df = pd.read_csv(in_file, sep='\t')
    df['mixing_id'] = df['mixing_ratio'].map(MR_MAP)

    if out_file:
        out_raw, f_fmt = os.path.splitext(out_file)
    else:
        out_raw, f_fmt = os.path.splitext(in_file)
    if f_fmt not in ('.pdf', '.png', '.jpg', '.jpeg', '.svg'):
        f_fmt = '.png'

    acc_MR_out = f'{out_raw}{f_fmt}'
    plot_accuracy_MR(df, acc_MR_out)

    # acc_MR_alg_out = f'{out_raw}_algorithms{f_fmt}'
    # plot_accuracy_MR_alg(df, acc_MR_alg_out)


def rename_mixing(x):
    p = [f'{float(i) * 100:.0f}' for i in x.split('-')]
    return ':'.join(p) + '%'


def plot_accuracy_MR(df, out_file):
    row_no = 1
    col_id = 'mixing_id'
    col_vals = df[col_id].unique()
    col_no = col_vals.size // row_no + col_vals.size % row_no + 1
    hue_id = 'algorithm'

    df['mixing_ratio_str'] = df['mixing_ratio'].apply(rename_mixing)
    df['doublet_rate_str'] = df['doublet_rate'].apply(lambda x: f'{x*100:.0f}%')

    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no * 15, row_no * 18))
    axes = np.reshape(axes, (row_no, col_no))
    for row in range(row_no):
        for col in range(col_no - 1): # Exclude legend col
            idx = row * col_no + col
            col_val = col_vals[idx]

            df_plot = df[df[col_id] == col_val]
            ax = axes[row][col]

            bp = sns.boxplot(
                data=df_plot,
                x='doublet_rate_str',
                y=METRIC,
                hue=hue_id,
                hue_order=sorted(df_plot[hue_id].unique()),
                palette=COLORS,
                ax=ax,
                fliersize=2,
                showfliers=False,
                linewidth=2,
                gap=0.1)
            sns.stripplot(
                data=df_plot,
                x='doublet_rate_str',
                y=METRIC,
                hue=hue_id,
                hue_order=sorted(df_plot[hue_id].unique()),
                palette=COLORS,
                ax=ax,
                linewidth=1,
                jitter=0.15,
                alpha=.8,
                size=6,
                dodge=True)

            ax.set_ylim((0.45, 1))
            ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1])

            # Plot linebreak
            # d = 0.02
            # y_br = 0.05
            # br_width = 0.02
            # ax.plot((-d, +d), (y_br-d, y_br+d),
            #     c='k', clip_on=False, transform=ax.transAxes, lw=2)
            # ax.plot((-d, +d), (y_br-d+br_width, y_br+d+br_width),
            #     c='k', clip_on=False, transform=ax.transAxes, lw=2)

            if col == 0:
                ax.set_ylabel('V-measure')
            else:
                ax.set_ylabel(None)
            ax.set_xlabel('Doublet rate')

            # Safe legend handles before removing legend
            handles, labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()

            y_label_pos = 1.05
            ax.set_title(f'{MR_MAP_reverse[col_val]}', y=y_label_pos)
        
            trans = mtransforms.ScaledTranslation(-1, 0, fig.dpi_scale_trans)
            ax.text(0, y_label_pos, PANELS[idx], transform=ax.transAxes + trans,
                fontsize='x-large', va='bottom', fontweight='bold', ha='center')

        # Add legend per row
        n_labels =  df_plot[hue_id].unique().size
        if hue_id == 'mixing_id':
            labels = [MR_MAP_reverse[int(i)] for i in labels[:n_labels]]
        if hue_id == 'algorithm':
            labels = [ALG_MAP[i] for i in labels[:n_labels]]
        else:
            labels = [format_MR(i) for i in labels[:n_labels]]

        axes[0][-1].axis('off')
        axes[0][-1].legend(
            handles[:n_labels], labels,
            title='Method',
            loc='center left',
            bbox_to_anchor=(-0.2, 0.5),
            ncol=1,
            borderaxespad=0.
        )

    plt.subplots_adjust(**MARGINS)
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


def format_MR(x):
    mr = [float(i) for i in x.split('-')]
    return ':'.join([f'{i*100:.0f}' for i in mr])


def plot_accuracy_MR_alg(df, out_file):
    row_id = 'algorithm'
    row_vals = df[row_id].unique()
    row_no = row_vals.size

    col_id = 'mixing_id'
    col_vals = df[col_id].unique()
    col_no = col_vals.size

    hue_id = 'mixing_ratio'

    df['mixing_ratio_str'] = df['mixing_ratio'].apply(rename_mixing)
    df['doublet_rate_str'] = df['doublet_rate'].apply(lambda x: f'{x*100:.0f}%')

    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no * 12, row_no * 12))
    axes = np.reshape(axes, (row_no, col_no))
    for row, row_val in enumerate(row_vals):
        for col, col_val in enumerate(col_vals):
            df_plot = df[(df[col_id] == col_val) & (df[row_id] == row_val)]
            ax = axes[row][col]

            bp = sns.boxplot(
                data=df_plot,
                x='doublet_rate_str',
                y=METRIC,
                hue=hue_id,
                hue_order=sorted(df_plot[hue_id].unique()),
                ax=ax,
                fliersize=2,
                showfliers=False,
                linewidth=1)
            sns.stripplot(
                data=df_plot,
                x='doublet_rate_str',
                y=METRIC,
                hue=hue_id,
                hue_order=sorted(df_plot[hue_id].unique()),
                ax=ax,
                linewidth=1,
                jitter=0.15,
                alpha=.8,
                size=6,
                dodge=True)

            ax.set_ylim((0, 1))
            ax.set_ylabel('V-measure')
            ax.set_xlabel('Doublet rate')

            handles, labels = ax.get_legend_handles_labels()
            n_labels =  df_plot[hue_id].unique().size
            if hue_id == 'mixing_id':
                labels = [MR_MAP_reverse[int(i)] for i in labels[:n_labels]]
            else:
                labels = [format_MR(i) for i in labels[:n_labels]]

            if n_labels == 1:
                n_col = 1
            else:
                n_col = 3

            ax.legend(handles[:n_labels], labels, title='Mixing ratio [%]',
                bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                ncol=n_col, mode="expand", borderaxespad=0.)

            idx = row * col_no + col
            trans = mtransforms.ScaledTranslation(-30/72, 1, fig.dpi_scale_trans)
            ax.text(0.0, 1.06, PANELS[idx], transform=ax.transAxes + trans,
                fontsize='xx-large', va='bottom', fontweight='bold', ha='center')

        ax2 = ax.twinx()
        ax2.set_ylabel(f'\nMethod:\n{ALG_MAP[row_val]}', fontsize='x-large',
            fontweight='bold')
        ax2.set_yticks([])
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        for tick in  ax.yaxis.majorTicks:
            tick.tick1line.set_markersize(0)

    fig.tight_layout()
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        help='Summary file from simulation runs.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for summary plot. default = <INPUT>.png')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output)