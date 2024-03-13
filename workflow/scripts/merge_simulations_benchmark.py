#!/usr/bin/env python3

import argparse
import os

import pandas as pd

def main(in_files, out_file):
    if not out_file:
        out_file = os.path.join(
            os.path.split(in_files[0])[0], 'benchmark_summary.tsv')
        
    if len(in_files) == 1 and os.path.isdir(in_files[0]):
        in_files = [os.path.join(in_files[0], i) for i in os.listdir(in_files[0])]

    df = pd.DataFrame()
    
    for in_file in in_files:
        if not in_file.endswith('txt'):
            continue
        file_split = os.path.split(in_file)[1].split('.')
        rep, alg = file_split[-3: -1]
        frac_str = '.'.join(file_split[:-3])
        frac, doublet = frac_str.split('-d')

        df_new = pd.read_csv(in_file, sep='\t')
        df_new['algorithm'] = alg
        df_new['mixing_ratio'] = frac
        df_new['doublet_rate'] = doublet
        df_new['replicate'] = rep

        df = pd.concat([df, df_new])

    df.to_csv(out_file, sep='\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Benchmark files from simulation runs.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for benchmarksummary plot. ' \
            'Default = <INPUT_DIR>/benchmark_summary.tsv')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output)