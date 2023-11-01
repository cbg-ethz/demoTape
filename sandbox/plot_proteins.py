#!/usr/bin/env python3

import argparse
import os

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_context('talk')


LABELS = {
    'HammingDist': 'Norm. Hamming Distance',
    'MP3': 'MP3 tree distance',
}
MARGINS = {
    'left': 0.1,
    'right': 0.95,
    'top': 0.95,
    'bottom': 0.05,
    'wspace': 0.4,
    'hspace': 0.6,
}
COL_NO = 5

NAME_MAP = {'0': '3867','1': '3841','2': '5962', '0+1': '3867+3841',
    '0+2': '3867+SID5962', '1+2': '3841+SID5962'}


def main(args):
    for i, input_file in enumerate(args.input):
        assign = pd.read_csv(args.assignment[i], index_col=0, dtype=str)
        cells = assign.columns.values.astype(np.bytes_)

        with h5py.File(input_file, 'r') as f:
            cell_ids = np.argwhere(
                np.isin(f['assays/protein_read_counts/ra/barcode'][:], cells)) \
                    .flatten()
            counts = f['assays/protein_read_counts/layers/read_counts'][cell_ids,:]
            AB_id = f['assays/protein_read_counts/ca/id'][()]
            # AB_id2 = f['assays/protein_read_counts/ca/antibody_id'][()]
            # AB_description = f['assays/protein_read_counts/ca/antibody_description'][()]

            df_wide = pd.DataFrame(counts, columns=AB_id.astype(str))

        df_wide['Barcode'] = assign.columns.values
        if len(args.input) > 1:
            import pdb; pdb.set_trace()
        else:
            df_wide['Cluster'] = assign.T['Cluster'].values

    try:
        df_wide['Cluster'] = df_wide['Cluster'].map(NAME_MAP)
    except:
        pass
    df_out = os.path.join(os.path.dirname(args.input[0]), 'protein_data.tsv')
    df_wide.to_csv(df_out, index=None,sep='\t')

    df_long = pd.melt(df_wide, id_vars=['Barcode', 'Cluster'])

    row_no = np.ceil(AB_id.size / 5).astype(int)

    fig, axes = plt.subplots(nrows=row_no, ncols=COL_NO,
        figsize=(row_no*6, COL_NO*10))

    for row in range(row_no):
        for col in range(COL_NO):
            AB_idx = row * COL_NO + col

            df_plot = df_long[df_long['variable'] == AB_id[AB_idx].astype(str)]
            ax = axes[row, col]
            bp = sns.boxplot(data=df_plot, x='variable', y='value', hue='Cluster',
                ax=ax, fliersize=2, showfliers=False, linewidth=1)
            sns.stripplot(data=df_plot, x='variable', y='value', hue='Cluster',
                ax=ax, linewidth=1, jitter=0.15, alpha=.2, size=4, dodge=True)

            ax.set_ylabel('log(Read Counts)')
            ax.set_xlabel('')
            ax.set_yscale('log')
            # if AB_idx != 0:
            #     bp.get_legend().remove()
            # else:
            handles, labels = ax.get_legend_handles_labels()
            new_len = int(len(handles) / 3)
            ax.legend(handles[:new_len], labels[:new_len], title='Cluster',
                bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                ncol=3, mode="expand", borderaxespad=0.)

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
        help='Input  h5 file containing protein information.')
    parser.add_argument('-a', '--assignment', type=str, nargs='+',
        help='Assignment file: cells -> cluster.')
    parser.add_argument('-o', '--output', type=str,
        help='Output file. If not set, print to stdout.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)