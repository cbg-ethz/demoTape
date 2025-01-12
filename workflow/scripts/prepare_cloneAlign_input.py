#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd


def main(args):
    if os.path.isfile(args.gene_annotation):
        ga_df = pd.read_csv(args.gene_annotation, sep='\t', index_col=0)
        # Drop duplicated genes: keep the ensemble ID with latest version number
        drop_idx = []
        for dup_gene in ga_df.index[ga_df.index.duplicated()]:
            dup_data = ga_df.loc[dup_gene]
            rm_id = np.argmin(
                dup_data['gene_id'].apply(lambda x: int(x.split('.')[1])).values)
            drop_idx.extend(
                np.argwhere((ga_df['gene_id'] == dup_data.iloc[rm_id, -1]).values)[0])
        ga_df = ga_df.iloc[~np.isin(np.arange(ga_df.shape[0]), drop_idx)]
        ga_df['gene_id'] = ga_df['gene_id'].apply(lambda x: x.split('.')[0])
        ga_df['id'] = ga_df['chrom'] + ga_df['arm']
        XY_genes = ga_df[(ga_df['chrom'] == 'chrX') | (ga_df['chrom'] == 'chrY')] \
            .index.values
        # Add gene alias for TCEB1
        if not 'TCEB1' in ga_df.index and 'ELOC' in ga_df.index:
            ga_df.loc['TCEB1'] = ga_df.loc['ELOC']
        if not 'FAM205A' in ga_df.index and 'SPATA31F1' in ga_df.index:
            ga_df.loc['FAM205A'] = ga_df.loc['SPATA31F1']
    else:
        raise IOError(f'Cannot open gene annotation: {args.gene_annotation}')

    CN_df = pd.read_csv(args.copy_numbers, index_col=0)
    # Drop  genes on sex chromosomes
    CN_df = CN_df.iloc[~np.in1d(CN_df.index, XY_genes)]
    # Drop genes without CNV
    CN_df = CN_df[(CN_df == 2).sum(axis=1) < 2]

    assert sum(~np.in1d(CN_df.index.values, ga_df.index.values)) == 0, \
        f'Not all gene names in {args.copy_numbers} are gencode names'
    
    if args.expression.endswith('.h5'):
        import h5py
        with h5py.File(args.expression, 'r') as f_exp:
            cells = f_exp['cell_attrs/cell_names'][:].astype(str)
            gene_ids = f_exp['gene_attrs/gene_ids'][:].astype(str)
            counts = f_exp['raw_counts'][:].astype(int)
    elif args.expression.endswith('.h5ad'):
        from anndata import read_h5ad
        adata = read_h5ad(args.expression)
        cells = adata.obs.index.values
        gene_ids = adata.var['gene_ids']
        counts = adata.X.toarray()
    exp_df = pd.DataFrame(counts, index=cells, columns=gene_ids)

    exp_df.columns = exp_df.columns.map(
        pd.Series(ga_df.index, index=ga_df['gene_id']).to_dict()).astype(str)
    # Drop  genes on sex chromosomes
    exp_df = exp_df.iloc[:, ~np.in1d(exp_df.columns, XY_genes)]

    exp_df = exp_df.iloc[:,exp_df.columns != 'nan']
    exp_genes = exp_df.columns

    if not args.output:
        out_base = os.path.dirname(args.copy_numbers)
        prefix = os.path.basename(args.copy_numbers).split('.')[0]
        args.output = os.path.join(out_base, prefix)
    elif args.output.endswith('.csv'):
        out_base = os.path.dirname(args.output)
        prefix = os.path.basename(args.output).split('.')[0]
        args.output = os.path.join(out_base, prefix)
    print('Writing files to: ' \
            f'{args.output}.<genes|chrArms>.cloneAlign_<expression|CNVclones>.csv')

    gene_overlap = sorted(
        CN_df.index.values[np.in1d(CN_df.index.values, exp_genes)])

    exp_gene_file = f'{args.output}.genes.cloneAlign_expression.csv'
    exp_df.loc[:, gene_overlap].T.to_csv(exp_gene_file)
    CN_gene_file = f'{args.output}.genes.cloneAlign_CNVclones.csv'
    CN_df.loc[gene_overlap].to_csv(CN_gene_file)


    print('The following genes were not found in the expression count matrix:')
    for i in CN_df.index[~np.in1d(CN_df.index.values, exp_genes)]:
        print(f'\t{i}')

    # Expand to all genes on same chromosome arm
    CN_df['id'] = CN_df.index.map(lambda x: ga_df.loc[x, 'id'])
    for chr_id, chr_id_data in CN_df.groupby('id'):
        chr_CN = int(round(chr_id_data['tumor'].mean()))
        if chr_CN == 2:
            continue
        chr_genes = ga_df[ga_df['id'] == chr_id].index.values
        chr_rel_genes = chr_genes[np.in1d(chr_genes, exp_genes)]

        new_exp = exp_df.loc[:, chr_rel_genes]
        new_CN = pd.DataFrame([[chr_CN, 2]] * chr_rel_genes.size,
            index=chr_rel_genes, columns=['tumor', 'healthy'])

        try:
            exp_df_chr = exp_df_chr.merge(new_exp, left_index=True, right_index=True)
            CN_df_chr = pd.concat([CN_df_chr, new_CN])
        except NameError:
            exp_df_chr = new_exp
            CN_df_chr = new_CN

    exp_gene_file = f'{args.output}.chrArms.cloneAlign_expression.csv'
    exp_df_chr.T.sort_index().to_csv(exp_gene_file)
    CN_gene_file = f'{args.output}.chrArms.cloneAlign_CNVclones.csv'
    CN_df_chr.sort_index().to_csv(CN_gene_file)
    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-exp', '--expression', type=str, required=True,
        help='Expression matrix (.h5).')
    parser.add_argument('-cn', '--copy_numbers', type=str, required=True,
        help='Copy number profiles (.csv).'),
    parser.add_argument('-ga', '--gene_annotation', type=str, required=True,
        default='resources/gencode.v19.annotation.cytoBand.tsv',
        help='Gencode gene-chromosome mapping (.tsv).')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output base for gene and cnv matrix files. ' \
            'Default = <CNV_FILE_DIR>/<CNV_FILE_BASE>')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)