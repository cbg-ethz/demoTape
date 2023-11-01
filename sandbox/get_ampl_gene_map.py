#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd

P2G_DEF = 'resources/Designer/3848-submitted.bed'
P2A_DEF = 'resources/Designer/3848-amplicon.bed'
OUT_DEF = 'resources/Designer/3848.amplicons2genes.bed'
SHM_DEF = 'resources/SHM_genes.csv'


def format_gene_name(x):
    if x.startswith('COSM'):
        return
    elif x.startswith('chr'):
        return
    else:
        return x.split(':')[0]


def main(args):
    df1 = pd.read_csv(args.position_to_gene, sep='\t', skiprows=1, header=None,
        names=['chr', 'start', 'end', 'gene_raw'])
    df1['gene'] = df1['gene_raw'].apply(format_gene_name)

    df1.set_index('chr', inplace=True)

    df2 = pd.read_csv(args.position_to_amplicon, sep='\t', header=None,
        names=['chr', 'start', 'end', 'amplicon'])
    df3 = pd.read_csv(args.shm_genes, true_values=['yes'], false_values=['no'], index_col=0)

    for i, data in df2.iterrows():
        chr_df = df1.loc[data["chr"],:]
        rel_pos = (chr_df['start'] >= data['start']) & (chr_df['end'] <= data['end'])

        if isinstance(rel_pos, (bool, np.bool_)):
            if rel_pos:
                df2.loc[i, 'gene'] = chr_df['gene']
                df2.loc[i, 'gene_raw'] = chr_df['gene_raw']
                if chr_df['gene'] in df3:
                    df2.loc[i, 'SHM'] = True
                else:
                    df2.loc[i, 'SHM'] = False
            else:
                continue
        else:
            pos_df = chr_df[rel_pos]

        if pos_df.size == 0:
            continue

        unique_genes = pos_df['gene'].value_counts()

        if unique_genes.size == 0:
            continue
        elif unique_genes.size == 1:
            df2.loc[i,'gene'] = pos_df.iloc[0,-1]
            df2.loc[i, 'gene_raw'] = pos_df.iloc[0,-2]
            if pos_df.iloc[0,-1] in df3:
                df2.loc[i, 'SHM'] = True
            else:
                df2.loc[i, 'SHM'] = False
        else:
            import pdb; pdb.set_trace()

    df2.to_csv(OUT_DEF, index=None)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p2g', '--position_to_gene', type=str, default=P2G_DEF,
        help='BED file containing: chr, pos1, pos2, gene.')
    parser.add_argument('-p2A', '--position_to_amplicon', type=str, default=P2A_DEF,
        help='BED file containing: chr, pos1, pos2, amplicon.')
    parser.add_argument('-shm', '--shm_genes', type=str, default=SHM_DEF,
        help='BED file containing SHM genes as first col.')
    parser.add_argument('-o', '--out_file', type=str, default=OUT_DEF,
        help='Output file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)