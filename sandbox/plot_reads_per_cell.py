#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main(args):
    df = pd.read_csv(args.doublets, sep='\t', index_col=0, usecols=[0, 2])
    reads = pd.read_csv(args.reads, sep=' ', header=None, index_col=1)
    reads.index += '-1'

    df['reads / cell'] = reads.loc[df.index]
    bp = sns.boxplot(data=df, x='doublet', y='reads / cell')
    # bp.set_yscale("log")
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--doublets', type=str, 
        help='COMPASS _cellAssignments.tsv file.')
    parser.add_argument('-r', '--reads', type=str, 
        help='missionbio .tube1.mapped.target.count.txt file.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)