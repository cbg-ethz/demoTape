#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
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
MR_MAP_reverse = {0: '33:33:33%', 1: '20:40:40%', 2: '20:50:30%', 3: '30:60:10%'}
ALG_MAP = {
    'scSplit': 'scSplit',
    'Euclidean': 'demoTape',
    'souporcell': 'souporcell',
    'demoTape': 'demoTape',
    'vireo': 'vireo'
}
PANELS = [f'{i}' for i in map(chr, range(65, 91))] # 97, 123 for small
COLORS = {
    'demoTape': '#F24C3D',
    'scSplit': '#F2BE22',
    'souporcell': '#22A699',
    'vireo': '#325FB0',
}


FONTSIZE = 30
DPI = 300
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'xtick.major.size': 2,
        'ytick.major.size': 2,
        'lines.linewidth': 2,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'large',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
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

    acc_MR_out = f'{out_raw}_accuracy_mixingRatio{f_fmt}'
    plot_accuracy_MR(df, acc_MR_out)

    acc_MR_alg_out = f'{out_raw}_accuracy_mixingRatio_algorithms{f_fmt}'
    plot_accuracy_MR_alg(df, acc_MR_alg_out)

    acc_out = f'{out_raw}_accuracy{f_fmt}'
    plot_accuracy(df, acc_out)


def plot_accuracy(df, out_file):
    hue_val = 'mixing_id'

    fig, ax = plt.subplots(figsize=(24, 12))
    bp = sns.boxplot(
        data=df,
        x='doublet_rate',
        y='accuracy',
        hue=hue_val,
        ax=ax,
        fliersize=2,
        showfliers=False,
        linewidth=1)
    sns.stripplot(
        data=df,
        x='doublet_rate',
        y='accuracy',
        hue=hue_val,
        ax=ax,
        linewidth=1,
        jitter=0.15,
        alpha=.8,
        size=6,
        dodge=True)

    ax.set_ylabel('Accuracy')
    ax.set_ylim((0, 1))
    ax.set_xlabel('Doublet rate')

    handles, labels = ax.get_legend_handles_labels()
    label_no = df[hue_val].nunique()
    if hue_val == 'mixing_id':
        labels = [MR_MAP_reverse[int(i)] for i in labels[:label_no]]
    else:
        labels = labels[:label_no]

    ax.legend(handles[:label_no], labels[:label_no], title='Mixing ratio',
        bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
        ncol=label_no // 2, mode="expand", borderaxespad=0.)
    plt.subplots_adjust(**{
        'left': 0.05,
        'right': 0.95,
        'top': 0.8,
        'bottom': 0.1,
        'wspace': 0.25,
        'hspace': 0.75,
    })
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


def rename_mixing(x):
    p = [f'{float(i) * 100:.0f}' for i in x.split('-')]
    return ':'.join(p) + '%'


def plot_accuracy_MR(df, out_file):
    row_no = 1
    col_id = 'mixing_id'
    col_vals = df[col_id].unique()
    col_no = col_vals.size // row_no + col_vals.size % row_no
    hue_id = 'algorithm'

    df['mixing_ratio_str'] = df['mixing_ratio'].apply(rename_mixing)
    df['doublet_rate_str'] = df['doublet_rate'].apply(lambda x: f'{x*100:.0f}%')

    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no * 12, row_no * 12))
    axes = np.reshape(axes, (row_no, col_no))
    for row in range(row_no):
        for col in range(col_no):
            idx = row * col_no + col
            col_val = col_vals[idx]

            df_plot = df[df[col_id] == col_val]
            ax = axes[row][col]

            bp = sns.boxplot(
                data=df_plot,
                x='doublet_rate_str',
                y='accuracy',
                hue=hue_id,
                ax=ax,
                fliersize=2,
                showfliers=False,
                linewidth=1)
            sns.stripplot(
                data=df_plot,
                x='doublet_rate_str',
                y='accuracy',
                hue=hue_id,
                ax=ax,
                linewidth=1,
                jitter=0.15,
                alpha=.8,
                size=6,
                dodge=True)

            # ax.set_ylim((0.3, 1))
            # ax.set_yticklabels([0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

            # # Plot linebreak
            # d = 0.02
            # y_br = 0.06
            # br_width = 0.02
            # ax.plot((-d, +d), (y_br-d, y_br+d),
            #     c='k', clip_on=False, transform=ax.transAxes, lw=2)
            # ax.plot((-d, +d), (y_br-d+br_width, y_br+d+br_width),
            #     c='k', clip_on=False, transform=ax.transAxes, lw=2)

            ax.set_ylabel('Accuracy')
            ax.set_xlabel('Doublet rate')

            handles, labels = ax.get_legend_handles_labels()
            label_no = df_plot[hue_id].nunique()
            if hue_id == 'mixing_id':
                labels = [MR_MAP_reverse[int(i)] for i in labels[:label_no]]
            if hue_id == 'algorithm':
                labels = [ALG_MAP[i] for i in labels[:label_no]]
            else:
                labels = [i.replace('-', ':') for i in labels[:label_no]]

            ax.set_title(f'Mixing ratio: {MR_MAP_reverse[col_val]}',
                fontweight='bold', y=1.3)

            if label_no == 6:
                ncol = 3
            else:
                ncol = label_no

            ax.legend(handles[:label_no], labels[:label_no], title='Method',
                bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                ncol=ncol // 2, mode="expand", borderaxespad=0.)

            trans = mtransforms.ScaledTranslation(-20/50, 1.5, fig.dpi_scale_trans)
            ax.text(0.0, 1.0, PANELS[idx], transform=ax.transAxes + trans,
                fontsize='xx-large', va='bottom', fontweight='bold', ha='center')

    fig.tight_layout()
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


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
                y='accuracy',
                hue=hue_id,
                ax=ax,
                fliersize=2,
                showfliers=False,
                linewidth=1)
            sns.stripplot(
                data=df_plot,
                x='doublet_rate_str',
                y='accuracy',
                hue=hue_id,
                ax=ax,
                linewidth=1,
                jitter=0.15,
                alpha=.8,
                size=6,
                dodge=True)

            ax.set_ylim((0, 1))
            ax.set_ylabel('Accuracy')
            ax.set_xlabel('Doublet rate')

            handles, labels = ax.get_legend_handles_labels()
            label_no = df_plot[hue_id].nunique()
            if hue_id == 'mixing_id':
                labels = [MR_MAP_reverse[int(i)] for i in labels[:label_no]]
            else:
                labels = [i.replace('-', ':') for i in labels[:label_no]]

            if label_no == 6:
                ncol = 3
            else:
                ncol = label_no

            ax.legend(handles[:label_no], labels[:label_no], title='Mixing ratio',
                bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                ncol=ncol, mode="expand", borderaxespad=0.)

            idx = row * col_no + col
            trans = mtransforms.ScaledTranslation(
                -20/72, 1.02 + 0.1 * (label_no // ncol), fig.dpi_scale_trans)
            ax.text(0.0, 1.0, PANELS[idx], transform=ax.transAxes + trans,
                fontsize='large', va='bottom', fontweight='bold')

        ax2 = ax.twinx()
        ax2.set_ylabel(f'\n\nMethod:\n{ALG_MAP[row_val]}', fontsize='x-large',
            fontweight='bold')
        ax2.set_yticks([])
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        for tick in  ax.yaxis.majorTicks:
            tick.tick1line.set_markersize(0)

    plt.subplots_adjust(**{
        'left': 0.05,
        'right': 0.9,
        'top': 0.9,
        'bottom': 0.1,
        'wspace': 0.25,
        'hspace': 0.33,
    })
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()



def plot_similarity(df, out_file):
    cmap='OrRd'

    col_no = 1
    row_no = 1
    fig, ax = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no*12, row_no*12))

    sns.set(font_scale=2.5)
    font_size = 30

    df_plot = pd.DataFrame()
    for name, val in df.items():
        t, p = name.split('->')
        if p in ['0+0', '1+1', '2+2', '3+3', '4+4', '5+5']:
            try:
                df_plot.loc[t, p[0]] += val
            except:
                df_plot.loc[t, p[0]] = val
        else:
            df_plot.loc[t, p] = val
    df_plot.fillna(0, inplace=True)
    df_plot = df_plot / df_plot.sum()

    sorted_idx = sorted(df_plot.index, key=lambda x: (x.count('+'), x[0]))
    sorted_cols = sorted(df_plot.columns, key=lambda x: (x.count('+'), x[0]))
    df_plot = df_plot.loc[sorted_idx, sorted_cols]

    hm = sns.heatmap(
        df_plot,
        annot=True,
        annot_kws={'size': font_size - 10},
        fmt='.3f',
        square=True,
        linecolor='lightgray',
        cmap=cmap,
        cbar_kws={'shrink': 0.5},
        linewidths=0.5,
        mask=df_plot == 0,
        ax=ax
    )

    ax.annotate(f'{dist} distance', xy=(0.5, 1.1), xytext=(0, 5), xycoords='axes fraction',
        textcoords='offset points', ha='center', va='baseline',
        annotation_clip=False)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=font_size)
    ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=font_size)
    ax.set_xlabel('Predicted', fontsize=font_size+10)
    ax.set_ylabel('True', fontsize=font_size + 10)

    MARGINS = {
        'left': 0.05,
        'right': 0.95,
        'top': 0.9,
        'bottom': 0.1,
        'wspace': 0.2,
        'hspace': 0.2,
    }
    plt.subplots_adjust(**MARGINS)

    if out_file:
        if not out_file.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg', '.svg')):
            out_file += '.png'
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        help='Summary file from simulation runs.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for summary tsv. default = <DIR:INPUT[0]>/summary.tsv')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output)