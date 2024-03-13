#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import tkinter as tk

from scipy.cluster.hierarchy import linkage, cut_tree
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context('talk')
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.transforms as mtransforms


ALG_MAP = {
    'scSplit': 'scSplit',
    'Euclidean': 'demoTape',
    'souporcell': 'souporcell',
    'demoTape': 'demoTape',
    'vireo': 'vireo'
}
COLORS = {
    'demoTape': '#F24C3D',
    'scSplit': '#F2BE22',
    'souporcell': '#22A699',
    'vireo': '#325FB0',
}
METRIC = 'cpu_time'


FONTSIZE = 30
DPI = 300
sns.set_style('whitegrid')
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
    if not out_file:
        out_file = os.path.join(os.path.split(in_file)[0], 'benchmark_summary.png')
    plot_runtime(df, out_file)


def plot_runtime(df, out_file):
    x_id = 'algorithm'

    fig, ax = plt.subplots(figsize=(12, 12))
    
    try:
        bp = sns.boxplot(
            data=df,
            x=x_id,
            order=sorted(COLORS.keys()),
            y=METRIC,
            ax=ax,
            palette=COLORS,
            fliersize=2,
            showfliers=False,
            linewidth=1)
    except tk.TclError:
        import matplotlib
        matplotlib.use('Agg')
        bp = sns.boxplot(
            data=df,
            x=x_id,
            order=sorted(COLORS.keys()),
            y=METRIC,
            ax=ax,
            palette=COLORS,
            fliersize=2,
            showfliers=False,
            linewidth=1)

    sns.stripplot(
        data=df,
        x=x_id,
        order=sorted(COLORS.keys()),
        y=METRIC,
        ax=ax,
        palette=COLORS,
        linewidth=1,
        jitter=0.15,
        alpha=.8,
        size=6,
        dodge=True)

    ax.set_ylim([-5, 200])
    ax.set_ylabel('CPU time [sec]')
    ax.set_xlabel('Algorithm')

    fig.tight_layout()
    if out_file:
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, 
        help='Benchmark summary files (.tsv) from simulation runs.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for benchmarksummary plot. ' \
            'Default = <INPUT_DIR>/benchmark_summary.png')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output)