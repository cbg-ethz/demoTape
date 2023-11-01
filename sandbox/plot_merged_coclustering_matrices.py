#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


SID2PAPER = {'SID3841': 'S1', 'SID3867': 'S2', 'SID5962': 'S3'}
FILE_EXT = 'png'
FIG_SIZE = (12, 5)
DPI = 300
PANELS = [f'{i})' for i in map(chr, range(97, 123))]
MARGINS = {
    'left': 0.05,
    'right': 0.99,
    'top': 0.9,
    'bottom': 0.05,
    'hspace': 0.2,
    'wspace': 0.2
}


def main(in_files, out_file):
    pngs = []
    for in_file in in_files:
        name_full = os.path.splitext(os.path.basename(in_file))[0]
        if name_full == 'sampling':
            continue
        sample = name_full.split('_')[0]
        if sample in SID2PAPER:
            sample = SID2PAPER[sample]
        alg = '_'.join(name_full.split('_')[3:])

        pngs.append([sample, alg, in_file])
    df = pd.DataFrame(pngs, columns=['sample', 'algorithm', 'filename'])

    row_id = 'sample'
    row_vals = df[row_id].unique()
    row_no = row_vals.size

    col_id = 'algorithm'
    col_vals = df[col_id].unique()
    col_no = col_vals.size

    fig, axes = plt.subplots(nrows=row_no, ncols=col_no,
        figsize=(FIG_SIZE[0]*col_no, FIG_SIZE[0]*row_no))
    axes = np.reshape(axes, (row_no, col_no))
    for row, row_val in enumerate(row_vals):
        for col, col_val in enumerate(col_vals):
            ax = axes[row, col]

            f_path = df[(df[row_id] == row_val) & (df[col_id] == col_val)]['filename'].iloc[0]
            arr_img = plt.imread(f_path,format='png')
            im = OffsetImage(arr_img)
            ab = AnnotationBbox(im, (1, 0), xycoords='axes fraction',
                box_alignment=(1,0), frameon=False, pad=0)
            ax.add_artist(ab)

    plt.subplots_adjust(**MARGINS)
    if not out_file:
        out_file = os.path.join(
            os.path.dirname(args.input[0]), f'coclustering_matrices_all.{FILE_EXT}')
        print(f'Writing file to: {out_file}')
    fig.savefig(out_file, dpi=DPI)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input directory (output dir of sampling run).')
    parser.add_argument('-o', '--out_file', type=str, default='',
        help='Output file. Default = <DIR: INPUT[0]>/coclustering_matrices_all.png')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.out_file)