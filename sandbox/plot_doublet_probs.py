#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_context('talk')
MARGINS = {
    'left': 0.15,
    'right': 0.9,
    'top': 0.95,
    'bottom': 0.2,
    'wspace': 0.25,
    'hspace': 0.75,
}

HIST_DEFAULT = {
    'alpha': 0.75,
    'fill': True,
    'binwidth': 0.05,
    'binrange': (0, 1),
    'element': 'bars',
    'stat': 'probability',
    'kde': False,
    'common_norm': False,
    'fill': True,
    'multiple': 'layer',
}


def main(args):
    df = pd.read_csv(args.assingment, sep='\t')
    dp = sns.histplot(df, x='doublet_probability', hue='doublet', **HIST_DEFAULT)

    dp.axes.set_xlabel('Doublet Probability')
    dp.axes.set_ylabel('Cell [%]')

    plt.subplots_adjust(**MARGINS)
    if not args.output:
        plt.show()
    else:
        if not args.output.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            args.output += '.png'
        dp.figure.savefig(args.output, dpi=300)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assingment', type=str,
        help='COMPASS _cellAssignments.tsv file.')
    parser.add_argument('-o', '--output', type=str, help='output file.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)