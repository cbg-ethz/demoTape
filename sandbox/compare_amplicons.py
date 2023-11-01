#!/usr/bin/env python3

import argparse
import json
import os
import re

import numpy as np
import pandas as pd

  

def main(args):
    for var_file in sorted(args.input_variants):
        df = pd.read_csv(var_file, dtype=str)
        df['idx'] = df['CHR'] + '_' + df['POS'] + '_' + df['REF'] + '_' + df['ALT']
        df.set_index('idx', inplace=True)

        called_amplicons = set(df['REGION'].values)

        for m_file in sorted(args.input_metrics):
            with open(m_file) as f:
                m = json.load(f)['additional_metrics']

            ampl_ids = ['high_expression_amplicons', 'failed_amplicons',
                'imbalanced_amplicons', 'low_performing_amplicons', 'noisy_amplicons']

            print(f'Files: {os.path.basename(var_file)} - {os.path.basename(m_file)} ')
            for i in ampl_ids:
                if i == 'low_performing_amplicons':
                    run_ampl = [i[0] for i in m[i]]
                else:
                    run_ampl = m[i]
                ol = called_amplicons.intersection(run_ampl)
                SNPs = df['REGION'].isin(list(ol)).sum()
                print(f'\toverlap {i}: {len(ol)} (#SNPs: {SNPs})')
        print()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-iv', '--input_variants', type=str, nargs='+',
        help='variant files to compare.')
    parser.add_argument('-im', '--input_metrics', type=str, nargs='+',
        help='Metric file to compare.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)