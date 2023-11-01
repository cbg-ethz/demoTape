#!/usr/bin/env python3

import argparse
import gzip
import numpy as np
import pandas as pd
import seaborn as sns


gt_map = {
    '0/0': 0, '0|0': 0, './.': 0, '.|.': 0,
    '0/1': 1, '0|1': 1, '1/0': 1, '1|0': 1,
    '0/2': 1, '0|2': 1, '2/0': 1, '2|0': 1,
    '0/3': 1, '0|3': 1, '3/0': 1, '3|0': 1,
    '1/1': 2, '1|1': 2, '2/2': 2, '2|2': 2, '3/3': 2, '3|3': 2,
    '1/2': 2, '1|2': 2, '2/1': 2, '2|1': 2,
    '1/3': 2, '1|3': 2, '3/1': 2, '3|1': 2,
    '2/3': 2, '2|3': 2, '3/2': 2, '3|2': 2,
    '3/4': 2, '3|4': 2, '4/3': 2, '4|3': 2,
}


def main(args):
    df = pd.read_csv(args.input, sep='\t', comment='#', header=None)

    line = '##'
    if args.input.endswith('.gz'):
        import gzip
        with gzip.open(args.input, 'r') as f:
            while True:
                line = f.readline().decode()
                if not line.startswith('##'):
                    break
    else:
        with open(args.input, 'r') as f:
            while True:
                line = f.readline()
                if not line.startswith('##'):
                    break
    df.columns = [i.split('-')[-1] for i in line.lstrip('#').rstrip().split('\t')]
    df.set_index(['CHROM', 'POS'], inplace=True)
    gt = df.iloc[:,7:].applymap(lambda x: gt_map[x.split(':')[0]])
    gt['unique'] = (gt != 0).sum(axis=1) == 1
    uniques = []
    for sample in sorted(df.columns[7:]):
        unique = gt[(gt[sample] > 0) & (gt['unique'])]
        uniques.append(unique)
        print(f'{sample}: {unique.shape[0]: >2}')
    gt_unique = pd.concat(uniques)

    plot(gt.iloc[:,:-1])
    plot(gt_unique.iloc[:,:-1])


def plot(data):
    cmap = sns.mpl.pyplot.get_cmap('inferno_r', 3)
    height = max(1, int(data.shape[0] // 5))
    width = max(1, int(data.shape[1] // 7))
    fig, ax = sns.mpl.pyplot.subplots(figsize=(width, height))
    sns.heatmap(data, linewidth=0.01, linecolor='black',
        square=True, cmap=cmap, vmin=0, vmax=2, ax=ax)
    ax.xaxis.tick_top()


    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)

    ax.set_yticks(np.arange(0.5, data.shape[0], 1))
    ax.set_xticks(np.arange(0.5, data.shape[1], 1))

    ax.set_xticklabels(data.columns, rotation=90, fontsize=8)
    y_labels = [f'{i[0]: >2}:{i[1]: >9}' for i in data.index.values]
    ax.set_yticklabels(y_labels, fontsize=8)

    sns.mpl.pyplot.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input vcf file.')
    parser.add_argument('-o', '--out_file', type=str, help='Output file.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)