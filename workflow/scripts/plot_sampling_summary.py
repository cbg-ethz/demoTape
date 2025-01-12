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
# Errorbars on barplots not defined for older versions
if [int(i) for i in sns.__version__.split('.')][1] < 12:
    raise ImportError(f'Requires seaborn version >= 0.12.X')

DEFAULT_CELLS = [500, 1000, 2000, 4000, 8000]
COLORS = {
    'BnpC': '#C25757', # red
    'SCG ( 8% dbt)': '#F0CF7F', # light cyan
    'SCG (25% dbt)': '#E9C61D', # cyan
    'COMPASS ( 8% dbt)': '#95AAD3', # light blue
    'COMPASS (25% dbt)': '#3A68AE', # blue
    '': '#FFFFFF'
}
ALG_ORDER = ['SCG ( 8% dbt)', 'SCG (25% dbt)', 'BnpC', 'COMPASS ( 8% dbt)',
    'COMPASS (25% dbt)']
COLORS_SNPS = {
    'SNPs_overlap': '#FFFFFF', # white
    'SNPs_missed': '#464546', # dark gray
    'SNPs_new': '#BBBBBC' # gray
}
MARGINS = {
    'left': 0.1,
    'bottom': 0.1,
    'right': 0.95,
    'top': 0.95,
    'wspace': 0.3,
    'hspace': 0.3 
}

SID2PAPER = {'SID3841': 'S1', 'SID3867': 'S2', 'SID5962': 'S3'}
MAX_CELLS = {'S1': 3732, 'S2': 7911, 'S3': 3677}
PANELS = [f'{i}' for i in map(chr, range(65, 91))] # 97, 123 for small


PANEL_SIZE = 10
FILE_EXT = 'png'
FONTSIZE = 30
DPI = 300
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'lines.linewidth': 2,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'xx-large',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 50
})
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['xtick.bottom'] = True


def get_x_ticks(df):
    cells = df['cells'].unique()
    pos = range(df['cells'].nunique())
    labels = [str(j) if i % 2 != 0 else '' for i, j in enumerate(cells)]
    return pos, labels


def add_panel_id(ax, label):
    trans = mtransforms.ScaledTranslation(-72/72, 0.0, ax.get_figure().dpi_scale_trans)
    ax.text(0, 1.0, label, transform=ax.transAxes + trans,
        fontsize='x-large', va='bottom', fontweight='bold', ha='center')


def add_alg_legend(ax):
    labels = ALG_ORDER
    ax.legend(
        [mpatches.Patch(color=COLORS[i], ec='black') for i in labels],
        labels,
        title='Algorithm',
        loc='center left',
        bbox_to_anchor=(-0.2, 0.5),
        ncol=1,
        borderaxespad=0.
    )


def add_snp_legend(ax):
    labels = ['both', 'only complete', 'only downsampled  ']
    colors = [COLORS_SNPS['SNPs_overlap'], COLORS_SNPS['SNPs_missed'],
        COLORS_SNPS['SNPs_new']]
    ax.legend(
        [mpatches.Patch(color=i, ec='black') for i in colors],
        labels,
        title='# SNPs detected in',
        loc='center left',
        bbox_to_anchor=(-0.2, 0.15),
        ncol=1,
        borderaxespad=0.
    )


# ----------------------------- first row --------------------------------------

def plot_SNP_no(df, ax, plot_xlabel=True):
    val_cols = ['SNPs_overlap', 'SNPs_missed', 'SNPs_new']
    id_cols = ['cells', 'run']
    rel_df = df[id_cols + val_cols].drop_duplicates()
   
    df_plot = pd.melt(rel_df, id_vars=id_cols, value_vars=val_cols)
    
    # Big green bar
    large_SNP_type = 'SNPs_overlap'
    bp = sns.barplot(
        data=df_plot[df_plot['variable'] == large_SNP_type],
        x='cells',
        y='value',
        errorbar='sd',
        capsize=0.2,
        lw=3,
        err_kws={'lw': 3},
        ec='black',
        ax=ax,
        color=COLORS_SNPS[large_SNP_type],
    )
    # Smaller red and yellow bars
    small_SNP_types = ['SNPs_missed', 'SNPs_new']
    bp = sns.barplot(
        data=df_plot[df_plot['variable'].isin(small_SNP_types)],
        x='cells',
        y='value',
        hue='variable',
        errorbar='sd',
        capsize=0.2,
        width = 0.9,
        gap=0.2,
        lw=3,
        ec='black',
        err_kws={'lw': 3},
        alpha=.75,
        ax=ax,
        palette=COLORS_SNPS,
    )

    if df_plot['cells'].nunique() > 6:
        x_ticks = get_x_ticks(df)
        ax.set_xticks(*x_ticks)
    if plot_xlabel:
        ax.set_xlabel('Sample size [cells]')
    else:
        ax.set_xlabel(None)

    y_max = df_plot[df_plot['variable'] == large_SNP_type]['value'].max()
    if y_max > 20:
        y_max_show = 56
    else:
        y_max_show = 21
    y_lim = (0, y_max_show)
    y_ticks = np.arange(y_lim[0], y_lim[1], 10)
    ax.set_yticks(y_ticks, np.abs(y_ticks))
    ax.set_ylim(y_lim)
    ax.set_ylabel('# SNPs')

    ax.get_legend().remove()


# ----------------------------- second row -------------------------------------

def plot_cluster_no(df, ax, plot_xlabel=True):
    y = 'Cluster'
    bp = sns.boxplot(
        data=df,
        x='cells',
        y=y,
        hue='algorithm',
        hue_order=ALG_ORDER,
        ax=ax,
        fliersize=2,
        linewidth=3,
        showfliers=False,
        palette=COLORS
    )
    sns.stripplot(
        data=df,
        x='cells',
        y=y,
        hue='algorithm',
        hue_order=ALG_ORDER,
        ax=ax,
        linewidth=1,
        jitter=0.15,
        alpha=.75,
        size=10,
        dodge=True,
        palette=COLORS
    )
    y_max = max(10, df[y].max())

    y_range = range(0, y_max + 1, 2)
    ax.set_yticks(y_range)
    ax.set_ylim((0, np.max(y_range) + 1))
    ax.set_ylabel('# Clusters')

    if df['cells'].nunique() > 6:
        x_ticks = get_x_ticks(df)
        ax.set_xticks(*x_ticks)
    if plot_xlabel:
        ax.set_xlabel('Sample size [cells]')
    else:
        ax.set_xlabel(None)

    ax.get_legend().remove()
    

# ----------------------------- third row --------------------------------------

def plot_metric(df, ax, metric='ARI', plot_xlabel=True):
    n_alg = df['algorithm'].nunique()
    bp = sns.boxplot(
        data=df,
        x='cells',
        y=metric,
        hue='algorithm',
        hue_order=ALG_ORDER,
        ax=ax,
        fliersize=2,
        showfliers=False,
        linewidth=3,
        palette=COLORS
    )
    sns.stripplot(
        data=df,
        x='cells',
        y=metric,
        hue='algorithm',
        hue_order=ALG_ORDER,
        ax=ax,
        linewidth=1,
        jitter=0.15,
        alpha=.75,
        size=10,
        dodge=True,
        palette=COLORS
    )

    if df['cells'].nunique() > 6:
        x_ticks = get_x_ticks(df)
        ax.set_xticks(*x_ticks)
    if plot_xlabel:
        ax.set_xlabel('Sample size [cells]')
    else:
        ax.set_xlabel(None)

    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_ylim((-0.05, 1.05))

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
        ax.set_ylabel(r'ARI')
    elif metric == 'Mutual info':
        ax.set_ylabel(r'Adjusted Mutual Information')
    elif metric == 'Norm. Hamming distance':
        ax.set_ylabel(r'Hamming distance (normalized)')
    
    ax.get_legend().remove()


def create_fig(samples, dfs, out_file):
    plot_fcts = [(plot_SNP_no, {}), (plot_cluster_no, {}),
        (plot_metric, {'metric': 'ARI'}),
        # (plot_metric, {'metric': 'Norm. Hamming distance'}),
    ]
    col_no = len(plot_fcts) + 1
    row_no = len(samples)

    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(PANEL_SIZE*col_no, PANEL_SIZE*row_no / 2))
    axes = np.reshape(axes, (row_no, col_no))

    for row, sample in enumerate(samples):
        df = dfs[row]
        try:
            name = SID2PAPER[sample]
        except KeyError:
            name = sample
        else:
            max_map = {df['cells'].max(): MAX_CELLS[name]}
            df.loc[:, 'cells'] = df['cells'].replace(max_map)

        for col, (plot_fct, plot_kwargs) in enumerate(plot_fcts):
            ax = axes[row][col]
            if col == 0:
                ax.annotate(f'{name}\n', xy=(0, 0.5),
                    xytext=(-ax.yaxis.labelpad - 50, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='xx-large', ha='right', va='center',
                )

            plot_fct(df, ax, plot_xlabel=row == row_no - 1,**plot_kwargs)
            add_panel_id(ax, PANELS[row * (col_no - 1) + col])

    # Add legends
    for row, _ in enumerate(samples):
        axes[row][-1].axis('off')
    add_snp_legend(axes[0][-1])
    add_alg_legend(axes[1][-1])

    plt.subplots_adjust(**MARGINS)
    print(f'Writing file to: {out_file}')
    fig.savefig(out_file, dpi=DPI)


def main(in_files, out_dir):
    if not out_dir:
        out_dir = os.path.dirname(in_files[0])

    dfs_in =  [pd.read_csv(i, sep='\t').drop_duplicates() for i in in_files]
    samples = [os.path.basename(i).split('_sampling_')[0] for i in in_files]

    for df in dfs_in:
        # Format doublet rates to percentage
        df.loc[:,'algorithm'] = df.apply(lambda x:
            f'{x["algorithm"]} ({x["doublet rate"] * 100:>2.0f}% dbt)' \
                if x["doublet rate"] > 0 else x['algorithm'],
            axis=1)

        # Select cell numbers that should be plotted, discard rest
        df.drop(df[~df['cells'].isin(args.cells)].index, inplace=True)
  
    for subset in dfs_in[0]['subset'].unique():
        subset='relevant'
        for min_cells in dfs_in[0]['min. #cells per cluster'].unique():
            out_file = os.path.join(out_dir,
                f'sampling_summary.min{min_cells}perCl_{subset}.{FILE_EXT}')
            dfs = [df[(df['min. #cells per cluster'] == min_cells) \
                & (df['subset'] == subset)] for df in dfs_in]

            create_fig(samples, dfs, out_file)


def get_files(in_dir):
    files = []
    for file in sorted(os.listdir(in_dir)):
        if file.endswith('sampling_summary.tsv'):
            files.append(os.path.join(in_dir, file))
    return files


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True,
        help='Directory containing <PREFIX>_sampling_summary.tsv files.'),
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory for result pngs. ' \
            f'Default = <DIR:INPUT[0]>/sampling_summary.min<X>perCl.{FILE_EXT}')
    parser.add_argument('-c', '--cells', type=int, nargs='+',
        default=DEFAULT_CELLS, help=f'Cell numbers to plot. ' \
            f'Default = {DEFAULT_CELLS}.')
    return parser.parse_args()


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.log:
            import sys
            sys.stderr = open(snakemake.log[0], 'w')
        main(snakemake.input, snakemake.output[0])
    else:
        args = parse_args()
        main(get_files(args.input), args.out_dir)