#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd


CHR_ORDER = {str(i): i for i in range(1, 23)}
CHR_ORDER.update({'X':23, 'Y':24})


def load_variants(in_file):
    return pd.read_csv(in_file, index_col=0)


def load_depth(in_file, whitelist=np.array([])):
    df = pd.read_csv(in_file, sep='\t', header=0, index_col=0)
    df.sort_index(inplace=True)

    def calc_gini(x):
        total = 0
        for i, xi in enumerate(x[:-1], 1):
            total += np.sum(np.abs(xi - x[i:]))
        return total / (len(x)**2 * np.mean(x))

    # Noisy amplicons (tapestri)
    gini = df.apply(lambda i: calc_gini(i.values))
    gini_threshold = gini.mean() + (2 * gini.std())
    noisy = gini > gini_threshold

    # low performing amplicons (tapestri)
    threshold = 0.2 * df.mean().mean()
    low_amplicons = df.mean() < threshold
    
    # Remove amplicons that are covered in less than 20% of all cells
    df = df.loc[:, (df == 0).sum() / df.shape[0] <= 0.8]
    # Remove amplicons with <= 1 read / cell
    low_cov = (df.sum() / df.shape[0]) < 1

    mean = df.mean().mean()
    # normalize read counts per cell
    df = df.apply(lambda x: (x / x.sum()), axis=1)
    # drop bad amplicons
    keep = df.columns.to_series().isin(whitelist)
    df = df.loc[:,~(noisy | low_amplicons | low_cov) | keep]

    # remove outliers
    df.clip(lower=df.quantile(0.1), upper=df.quantile(0.9), axis=1, inplace=True)
    # normalize amplicons
    df = df / df.mean()

    return df.T, mean


def load_panel(in_file):
    df = pd.read_csv(in_file, comment='#', sep='\t', header=None, index_col=-1)
    df[0] = df[0].str.replace('chr', '')

    df['locus'] = df[0].astype(str) + ':' + df[1].astype(str) \
            + '-' + df[2].astype(str)
    no_gene = df[df[3] == '.'].index
    df.loc[no_gene, 3] = df.loc[no_gene, 'locus']
    return df


def load_gene_annoptation(in_file):
    df = pd.read_csv(args.gene_annotation, sep='\t', index_col=0)
    # Drop duplicated genes: keep the ensemble ID with latest version number
    drop_idx = []
    for dup_gene in df.index[df.index.duplicated()]:
        dup_data = df.loc[dup_gene]
        rm_id = np.argmin(
            dup_data['gene_id'].apply(lambda x: int(x.split('.')[1])).values)
        drop_idx.extend(
            np.argwhere((df['gene_id'] == dup_data.iloc[rm_id, -1]).values)[0])
    df = df.iloc[~np.isin(np.arange(df.shape[0]), drop_idx)]
    df['gene_id'] = df['gene_id'].apply(lambda x: x.split('.')[0])
    df['id'] = df['chrom'] + df['arm']

    return df


def load_cell_annotation(in_file):
    df = pd.read_csv(in_file, index_col=0)
    if 'assignmnet' in df.columns:
        df.rename({'assignmnet': 'assignment'}, axis=1, inplace=True)
    return df


def normalize_depth(df_d, df_ca):
    healthy_cl = df_ca[df_ca['assignment'] == 'healthy'].iloc[0]['cluster']
    healthy_cells = df_ca['cluster'] == healthy_cl

    # Normalize such that avg. healthy cells depth = 2
    return df_d.apply(lambda x: x / x[healthy_cells].mean() * 2, axis=1)


def save_data(df_v, df_d, out_pattern, region):
    var_file = out_pattern.format(region=region, dtype='variants')
    print(f'Writing file: {var_file}')
    df_v.to_csv(var_file)

    # Sort by chromosome
    idx_sort = sorted(df_d.index, key=lambda i: (CHR_ORDER[i.split('_')[0]], i))
    df_d = df_d.loc[idx_sort]

    depth_file = out_pattern.format(region=region, dtype='regions')
    print(f'Writing file: {depth_file}')
    df_d.round().astype(int).to_csv(depth_file, header=False)


def update_chrArm_map(df_ga, gene_chrArm, no_genes):
    for pos_id in no_genes:
        chrom, pos_start_end = pos_id.split(':')
        pos_mean = np.mean([int(i) for i in pos_start_end.split('-')])

        df_chr = df_ga[(df_ga['chrom'] == f'chr{chrom}')]
        overlap = df_chr[(pos_mean >= df_chr['start']) \
            & (pos_mean <= df_chr['end'])]
        if overlap.shape[0] == 0:
            min_pos =  df_chr[(pos_mean >= df_chr['start'])]
            max_pos = df_chr[(pos_mean <= df_chr['end'])]
            if min_pos['arm'].nunique() == 1:
                arm = min_pos['arm'].unique()[0]
            elif max_pos['arm'].nunique() == 1:
                arm = max_pos['arm'].unique()[0]
            else:
                arm = pos_id
                print(f'Cannot map position {pos_id} to chromosome arm')
        else:
            overlap_arm = overlap['arm'].unique()
            if len(overlap_arm) == 1:
                arm = overlap_arm[0]
            else:
                arm = pos_id
                print(f'Cannot map position {pos_id} to chromosome arm')
        gene_chrArm[pos_id] = f'chr{chrom}{arm}'


def main(args):
    if not args.output:
        args.output = args.variants.split('variants.csv')[0]
    
    if 'filtered_' in args.variants:
        filter_str = '.filtered'
        args.output = args.output.split('filtered_')[0].rstrip('.')
    else:
        filter_str = ''
    out_pattern = f'{args.output}.{{region}}{filter_str}_{{dtype}}.csv'

    df_v = load_variants(args.variants)
    df_d, avg_d = load_depth(args.depths, whitelist=df_v['REGION'].unique())
    # Make sure that cell order is the same
    df_d = df_d[df_v.columns[6:]]
    
    df_p = load_panel(args.panel)

    if args.cell_annotation:
        df_ca = load_cell_annotation(args.cell_annotation)
        # Remove doublets
        if 'doublets' in df_ca['assignment'].unique():
            dbt_cells = df_ca[df_ca['assignment'] == 'doublets'].index
            df_v.drop(dbt_cells, axis=1, inplace=True)
            df_d.drop(dbt_cells, axis=1, inplace=True)
        df_d = normalize_depth(df_d, df_ca)

    # change read depth to resemble original avg mean
    df_d = df_d * avg_d / df_d.mean().mean()

    df_d['CHR'] = df_d.index.map(dict(zip(df_p.index, df_p[0])))
    ampl_index = df_d.apply(lambda x: f'{x.CHR}_{x.name}', axis=1)

    # Safe amplicon data
    df_v['REGION'] = df_v['NAME']

    save_data(df_v, df_d.rename(ampl_index).drop('CHR', axis=1), out_pattern,
        'amplicons')
    
    # Aggregate over genes
    ampl_gene = dict(zip(df_p.index, df_p[3]))
    df_v['REGION'] = df_v['NAME'].map(ampl_gene)

    df_d.index = df_d.index.map(ampl_gene)
    df_d.index.name = 'gene'
    df_d = df_d.groupby('gene').sum()
    df_d['CHR'] = df_d.index.map(dict(zip(df_p[3], df_p[0])))
    gene_index = df_d.apply(lambda x: f'{x.CHR}_{x.name}', axis=1)
    save_data(df_v, df_d.rename(gene_index).drop('CHR', axis=1), out_pattern,
        'genes')

    if not args.gene_annotation:
        exit()

    # Aggregate over chromosome arms
    df_ga = load_gene_annoptation(args.gene_annotation)
    # Add gene alias for TCEB1
    if not 'TCEB1' in df_ga.index and 'ELOC' in df_ga.index:
        df_ga.loc['TCEB1'] = df_ga.loc['ELOC']
    if not 'FAM205A' in df_ga.index and 'SPATA31F1' in df_ga.index:
        df_ga.loc['FAM205A'] = df_ga.loc['SPATA31F1']

    # sanity check panel genes to be also in Gencode DB
    p_genes = df_p[3].unique()
    unique_genes = p_genes[np.argwhere(~np.isin(p_genes, df_ga.index)).ravel()]
    unique_genes = [i for i in unique_genes if not (':' in i and '-' in i)]
    print(f'Genes not found in Gencode DB (RNA genes?): {unique_genes}')

    gene_chrArm = dict(zip(df_ga.index, df_ga['id']))
    no_genes = df_d[df_d.index.str.contains(':')].index.values
    update_chrArm_map(df_ga, gene_chrArm, no_genes)

    df_v['REGION'] = df_v['REGION'].map(gene_chrArm)

    df_d.index = df_d.index.map(gene_chrArm)
    df_d.index.name = 'chrArm'
    df_d = df_d.groupby('chrArm').sum()
    df_d['CHR'] = [i[3:-1] for i in df_d.index]
    chrArm_index = df_d.apply(lambda x: f'{x.CHR}_{x.name}', axis=1)
    save_data(df_v, df_d.rename(chrArm_index).drop('CHR', axis=1), out_pattern,
        'chrArms')
        

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--variants', type=str, required=True,
        help='Filtered variant file (from mosaic_preprocessing.py) (.csv).'),
    parser.add_argument('-d', '--depths', type=str, required=True,
        help='Reads per barcode distribution from Tapestri processing (.tsv).'),
    parser.add_argument('-p', '--panel', type=str, required=True,
        help='Annotated Tapestri panel bed file')
    parser.add_argument('-ga', '--gene_annotation', type=str,
        help='Gencode gene-chromosome mapping (.tsv).')
    parser.add_argument('-ca', '--cell_annotation', type=str, default='',
        help='Cell order (barcodes) in csv format.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output base name.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)