#!/usr/bin/env python3

import argparse
import os
import shutil

import loompy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
import umap
sns.set_context('talk')


SAMPLE_MAP = {'SID3841': 'S1', 'SID3867': 'S2', 'SID5962': 'S3',
    'mixed_SID3867_SID5962_SID3841': 'MS1'}
Cluster2PAPER = {'0': 'S3', '1': 'S1', '2': 'S2',
    '0+1': 'S1+S3', '1+2': 'S1+S2', '0+2': 'S2+S3'}


def get_variants_data(in_files):
    for in_file in in_files:
        name = os.path.basename(in_file).split('.cells')[0]
        if name in SAMPLE_MAP:
            name = SAMPLE_MAP[name]

        df_new = pd.read_csv(in_file)
        df_new['idx'] = df_new['CHR'].astype(str) + '_' + df_new['POS'].astype(str) \
            + '_' + df_new['REF'] + '_' + df_new['ALT']
        df_new.set_index('idx', drop=True, inplace=True)
        df_new.drop(['CHR', 'POS', 'REF', 'ALT', 'REGION', 'NAME', 'FREQ'], axis=1, inplace=True)

        df_new = df_new.applymap(lambda x: float(x.split(':')[2]))

        df_new.replace(3, np.nan, inplace=True)
        if name != 'MS1':
            df_new.columns = df_new.columns + f'.{name}'
            df_new = df_new.T.fillna(df_new.T.mean()).T

        cell2sample_new = pd.Series(name, index=df_new.columns)

        try:
            df = df.merge(df_new, left_index=True, right_index=True, how='inner')
            cell2sample = cell2sample.append(cell2sample_new)
        except NameError:
            df = df_new
            cell2sample = cell2sample_new

    return df, cell2sample


def main(args):
    df, c2s_map = get_variants_data(args.input)
    if args.assignment:
        c2s_map_new = pd.read_csv(
            args.assignment, index_col=0, dtype='str', sep='\t').T.squeeze()
        if not c2s_map_new.map(Cluster2PAPER).isna().any():
            c2s_map_new = c2s_map_new.map(Cluster2PAPER)

        c2s_map.loc[c2s_map_new.index] = c2s_map_new
        for i in c2s_map_new.unique():
            cells = c2s_map_new[c2s_map_new == i].index
            cells = [i for i in cells if i in df.columns.values]
            df[cells] = df[cells].T.fillna(df[cells].T.mean()).T
    else:
        df = df.T.fillna(df.T.mean()).T

    samples = np.unique(c2s_map.values)
    try:
        samples = sorted(samples, key = lambda x: (len(x), int(x[0])))
    except ValueError:
        pass
    cell_colors = c2s_map.map(dict(zip(samples, range(len(samples)))))

    plot_df = df.values.T# np.nan_to_num(df.values.T, nan=1)
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(plot_df)

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.scatter(
        embedding[:, 0],
        embedding[:, 1],
        s=5,
        c=[sns.color_palette()[x] for x in cell_colors],
        alpha=0.5)
    ax.set_aspect('equal', 'datalim')

    ax.set_title('Multiplexed scDNA-seq data', fontsize=24)
    ax.set_ylabel('UMAP 1')
    ax.set_xlabel('UMAP 2')

    handles = [mpatches.Patch(color=sns.color_palette()[x]) for x in range(len(samples))]
    ax.legend(handles, samples, frameon=True, title='Sample')

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
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input variants.tsv file(s) of multiplexed samples.')
    parser.add_argument('-a', '--assignment', type=str, default='',
        help='Assignment file: cells -> cluster.')
    parser.add_argument('-o', '--output', type=str,
        help='Output file. If not set, print to stdout.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)