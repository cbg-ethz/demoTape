#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial import distance
import seaborn as sns
from tqdm import tqdm

sns.set_context('talk')


def merge_gt(gt_in):
    # 00=0, 11=1, 22=2, 33=3; 03=0, 13=1, 23=2; 01=1, 02=1, 12=1
    if gt_in[0] == gt_in[1]:
        return gt_in[0]
    elif np.isnan(gt_in[0]):
        return gt_in[1]
    elif np.isnan(gt_in[1]):
        return gt_in[0]
    else:
        return 1


def get_COMPASS_results(assign_file):
    cells_att = pd.read_csv(assign_file, sep='\t', index_col=0)

    gt_file = assign_file.replace('cellAssignments', 'nodes_genotypes')
    gt_nodes = pd.read_csv(gt_file, sep='\t')
    gt_nodes.set_index(gt_nodes['node'].str.replace('Node ', '').astype(int),
        inplace=True)
    gt_nodes.drop('node', axis=1, inplace=True)

    sgt_idx = cells_att[cells_att['doublet'] == 'no'].index
    sgt_assign = cells_att.loc[sgt_idx, 'node'].astype(int)
    gt = gt_nodes.loc[sgt_assign]
    gt.set_index(sgt_idx, inplace=True)

    dbt_idx = cells_att[cells_att['doublet'] == 'yes'].index
    dbt_cl_df = cells_att.loc[dbt_idx, 'node'].astype(str).str.split('+') \
        .apply(lambda x: [int(i) for i in x])
    dbt_assign = np.concatenate(dbt_cl_df.values).reshape((dbt_idx.size, 2))
    dbt_gt = pd.DataFrame(
        np.apply_along_axis(merge_gt, axis=1, arr=gt_nodes.values[dbt_assign]),
        index=dbt_idx, columns=gt_nodes.columns)

    assign_df = pd.DataFrame([sgt_assign.values + 1, np.zeros(sgt_assign.size)], columns=sgt_idx).T
    assign_df = assign_df.append(pd.DataFrame(dbt_assign + 1, index=dbt_idx))

    return assign_df.loc[cells_att.index], gt.append(dbt_gt).loc[cells_att.index]


def get_BnpC_results(assign_file, cl):
    assign = np.array(
        pd.read_csv(assign_file, sep='\t').loc[0, 'Assignment'].split(' '),
        dtype=int)

    gt_file = os.path.join(os.path.dirname(assign_file), 'genotypes_posterior_mean.tsv')
    gt = pd.read_csv(gt_file, sep='\t', index_col=0).T

    path_dirs = os.path.abspath(assign_file).split(os.sep)
    base_dir = os.path.join(*[os.sep] + path_dirs[:-4])

    sample_file = [os.path.join(base_dir, i) for i in os.listdir(base_dir)\
        if f'cells.{cl}_variants.csv' in i][0]
    with open(sample_file, 'r') as f:
        cells = f.readline().strip().split(',')[7:]
    gt.index = cells
    return pd.Series(assign, index=cells), gt


def get_SCG_results(assign_file):
    cells_att = pd.read_csv(assign_file, sep='\t', index_col=0)

    gt_file = os.path.join(os.path.dirname(assign_file), 'Cluster_genotypes.tsv')
    gt_nodes = pd.read_csv(gt_file, sep='\t', index_col=0)

    gt = gt_nodes.loc[cells_att.loc[:, 'cluster_id']]
    return cells_att['cluster_id'], gt.set_index(cells_att.index)


def plot_coclustering_matrix(df, out_file=''):
    Z = linkage(df, 'ward', 'euclidean')
    cell_order = leaves_list(Z)

    fig, ax = plt.subplots(figsize=(15, 15))
    hm = sns.heatmap(
        df.iloc[cell_order,cell_order],
        ax=ax,
        square=True,
        vmin=0,
        vmax=1,
        cmap='Reds',
        linewidths=0,
        linecolor='lightgray',
        cbar_kws={'label': 'avg. pairwise-cluster-distance'}
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel('Cells')
    ax.set_xlabel('Cells')

    plt.subplots_adjust(**{
        'left': 0.05,
        'right': 0.95,
        'top': 0.95,
        'bottom': 0.05
    })
    if out_file:
        if not out_file.lower().endswith(('.pdf', '.png', '.jpg', '.jpeg')):
            out_file += '.png'
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()


def get_data(in_files):
    data = {}
    for in_file in in_files:
        if re.search('COMPASS', in_file):
            alg = 'COMPASS'
            dbt = re.search('/r\d+\.d(0.\d+)_', in_file).group(1)
            cl = re.search('COMPASS/cl(\d+)/', in_file).group(1)
            assign, gt = get_COMPASS_results(in_file)
        elif re.search('BnpC', in_file):
            alg = 'BnpC'
            dbt = '-1'
            cl = re.search('BnpC/cl(\d+)/', in_file).group(1)
            assign, gt = get_BnpC_results(in_file, cl)
        elif re.search('SCG', in_file):
            alg = 'SCG'
            dbt = re.search('/r\d+\.d(0.\d+)/', in_file).group(1)
            cl = re.search('SCG/cl(\d+)/', in_file).group(1)
            assign, gt = get_SCG_results(in_file)

        if not cl in data:
            data[cl] = {}
        if not alg in data[cl]:
            data[cl][alg] = {}
        if not dbt in data[cl][alg]:
            data[cl][alg][dbt] = {'assign': [], 'geno': []}

        data[cl][alg][dbt]['assign'].append(assign)
        data[cl][alg][dbt]['geno'].append(gt)

    gr_truth = {}
    for cl, cl_data in data.items():
        gr_truth[cl] = {}
        for alg, alg_data in cl_data.items():
            gr_truth[cl][alg] = {}
            for dbt, dbt_data in alg_data.items():
                gr_truth[cl][alg][dbt] = {
                    'assign': get_coclustering_matrix(dbt_data['assign']),
                    'geno': get_geno_mean(dbt_data['geno'])}

    return gr_truth


def get_geno_mean(genos):
    rows = genos[0].index
    cols = genos[0].columns
    vals = []
    for geno in genos:
        vals.append(geno.loc[rows, cols].values)
    return pd.DataFrame(np.stack(vals).mean(axis=0), index=rows, columns=cols)


def get_coclustering_matrix(assigns):
    cells = assigns[0].index
    vals = []
    for assign_raw in assigns:
        assign = assign_raw.loc[cells].values
        if assign.ndim == 2:
            vals.append((assign[None, :,0] == assign[:, None, 0]).astype(int) \
                + (assign[None, :,0] == assign[:, None, 1]).astype(int))
        else:
            vals.append((assign[None,:] == assign[:,None]).astype(int))
    return pd.DataFrame(np.stack(vals).mean(axis=0), index=cells, columns=cells)



def main(in_files, out_dir):
    data = get_data(in_files)

    for cl, cl_data in data.items():
        for alg, alg_data in cl_data.items():
            for dbt, dbt_data in alg_data.items():
                out_file_cocl = os.path.join(out_dir,
                    f'cl{cl}_coclustering_{alg}_d{dbt}.png')
                plot_coclustering_matrix(dbt_data['assign'], out_file_cocl)



def get_COMPASS_files(dir_path):
    files = []
    for in_file in os.listdir(dir_path):
        if not in_file.endswith('_cellAssignments.tsv'):
            continue
        files.append(os.path.join(dir_path, in_file))
    return files


def get_BnpC_files(dir_path):
    files = []
    for run in os.listdir(dir_path):
        files.append(os.path.join(dir_path, run, 'assignment.txt'))
    return files


def get_SCG_files(dir_path):
    files = []
    for run in os.listdir(dir_path):
        files.append(os.path.join(dir_path, run, 'assignments.tsv'))
    return files


def get_files(in_dir):
    files = []
    for sub_dir in os.listdir(in_dir):
        dir_path = os.path.join(in_dir, sub_dir)
        if sub_dir == 'samples':
            continue
        elif sub_dir == 'COMPASS':
            files.extend(get_COMPASS_files(dir_path))
        elif sub_dir == 'BnpC':
            files.extend(get_BnpC_files(dir_path))
        elif sub_dir == 'SCG':
            files.extend(get_SCG_files(dir_path))
    return files


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        help='Input directory (output dir of sampling run).')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output file. Default = <INPUT>/cooccurence_matrix.png')
    return parser.parse_args()


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.log:
            import sys
            sys.stderr = open(snakemake.log[0], 'w')
        main(snakemake.input, snakemake.params.out_dir)
    else:
        args = parse_args()
        main(get_files(args.input), args.out_dir)