#!/usr/bin/env python3

import argparse
import os
import shutil
import tempfile

import loompy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_context('talk')


def get_assignment_dict(assignment_file):
    clusters = {}
    assign = pd.read_csv(assignment_file, index_col=0, dtype='str', sep='\t')
    for cl in np.unique(assign):
        clusters[cl] = assign.columns[(assign == cl).values.flatten()].values
    return clusters


def get_loom_data(in_file):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_in_file = os.path.join(temp_dir, os.path.basename(in_file))
        shutil.copy2(args.input, temp_in_file)
        with loompy.connect(temp_in_file) as ds:
            cells = ds.col_attrs['barcode'].astype(str)
            DP = ds.layers['DP'][()]
            GQ = ds.layers['GQ'][()]
            AD = ds.layers['AD'][()]
            RO = ds.layers['RO'][()]
    VAF = np.nan_to_num(AD / DP, nan=0)
    return {'DP': DP, 'GQ': GQ, 'AD': AD, 'RO': RO, 'VAF': VAF}, cells


def get_variants_data(in_file):
    df = pd.read_csv(in_file, index_col=[0, 1], dtype={'CHR': str})
    cells = df.columns[5:].values
    DP = df.iloc[:,5:].applymap(lambda x: float(x.split(':')[0])).values
    AD = df.iloc[:,5:].applymap(lambda x: float(x.split(':')[1])).values
    VAF = np.nan_to_num(AD / DP, nan=0)
    return {'DP': DP, 'AD': AD, 'VAF': VAF}, cells


def main(args):
    if args.input.endswith('.loom'):
        data, cells = get_loom_data(args.input)
    elif args.input.endswith('.csv'):
        data, cells = get_variants_data(args.input)

    assign = get_assignment_dict(args.assignment)

    row_no = 1
    col_no = len(data)

    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(col_no*12, row_no*12))

    for i, (data_str, data_col) in enumerate(data.items()):
        ax = axes[i]

        df_plot = []
        df_labels = []

        print(f'{data_str}:')
        for cl_id in ['0', '1', '2', '0+1', '0+2', '1+2']:
            cl_cells = assign[cl_id]
            cl_idx = np.argwhere(np.in1d(cells, cl_cells)).flatten()
            cl_data = data_col[:,cl_idx].flatten()
            df_plot.append(cl_data)
            df_labels.append(cl_id)
            print(f'\t{cl_id}:\t{cl_data.mean():.3f} +/- {cl_data.std():.3f}')

        bp = sns.boxplot(data=df_plot, ax=ax, fliersize=2, linewidth=1)

        ax.set_ylabel(data_str)
        ax.set_xticklabels(df_labels)
        ax.set_xlabel('Cluster')

    MARGINS = {
        'left': 0.1,
        'right': 0.95,
        'top': 0.95,
        'bottom': 0.1,
        'wspace': 0.5,
    }
    plt.subplots_adjust(**MARGINS)
    if not args.output:
        plt.show()
    else:
        if not args.output.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            args.output += '.png'
        fig.savefig(args.output, dpi=300)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        help='Input loom file of multiplexed samples.')
    parser.add_argument('-a', '--assignment', type=str,
        help='Assignment file: cells -> cluster.')
    parser.add_argument('-o', '--output', type=str,
        help='Output file. If not set, print to stdout.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)