#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

FRAC_PATTERN = re.compile(r'/([\d\.]+-[\d\.]+-[\d\.]+)')
DBT_PATTERN = re.compile(r'-d(0\.\d{2})')
DEPTH_PATTERN = re.compile(r'-dp(0\.\d{2})')


def main(in_files, out_file):
    for in_file in in_files:
        df_new = pd.read_csv(in_file, sep='\t')

        df_new['doublet_rate'] = float(DBT_PATTERN.search(in_file).group(1))
        df_new['mixing_ratio'] = FRAC_PATTERN.search(in_file).group(1)
        try:
            depth = DEPTH_PATTERN.search(in_file).group(1)
        except (KeyError, AttributeError):
            depth = 1
        df_new['depth_downsample'] = depth

        try:
            df = pd.concat([df, df_new], axis=0)
        except NameError:
            df = df_new

    df = df[['mixing_ratio', 'doublet_rate', 'depth_downsample', 'algorithm',
        'rep', 'v_measure']]

    if not out_file:
        out_file = os.path.join(os.path.dirname(in_files[0]), 'summary.tsv')

    df.to_csv(out_file, sep='\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Summary files from simulation runs.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for summary tsv. default = <DIR:INPUT[0]>/summary.tsv')
    return parser.parse_args()


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.log:
            import sys
            sys.stderr = open(snakemake.log[0], 'w')
        main(snakemake.input, snakemake.output[0])
    else:
        args = parse_args()
        main(args.input, args.output)