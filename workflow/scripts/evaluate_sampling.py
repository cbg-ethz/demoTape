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
from scipy.spatial.distance import pdist, euclidean
from sklearn.metrics.cluster import adjusted_rand_score, adjusted_mutual_info_score
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
    assign_df = pd.DataFrame([sgt_assign.values + 1, np.zeros(sgt_assign.size)],
        columns=sgt_idx).T

    dbt_idx = cells_att[cells_att['doublet'] == 'yes'].index
    if dbt_idx.size > 0:
        dbt_cl_df = cells_att.loc[dbt_idx, 'node'].astype(str).str.split('+') \
            .apply(lambda x: [int(i) for i in x])
        dbt_assign = np.concatenate(dbt_cl_df.values).reshape((dbt_idx.size, 2))
        dbt_gt = pd.DataFrame(
            np.apply_along_axis(merge_gt, axis=1, arr=gt_nodes.values[dbt_assign]),
            index=dbt_idx, columns=gt_nodes.columns)
        assign_df = assign_df.append(pd.DataFrame(dbt_assign + 1, index=dbt_idx))
        gt = gt.append(dbt_gt).loc[cells_att.index]

    return assign_df.loc[cells_att.index], gt


def get_BnpC_results(assign_file):
    assign = np.array(
        pd.read_csv(assign_file, sep='\t').loc[0, 'Assignment'].split(' '),
        dtype=int)

    gt_file = os.path.join(os.path.dirname(assign_file), 'genotypes_posterior_mean.tsv')
    gt = pd.read_csv(gt_file, sep='\t', index_col=0).T

    path_dirs = os.path.abspath(assign_file).split(os.sep)
    sample_file = os.path.join(*[os.sep] + path_dirs[:-3] \
        + ['samples', f'{path_dirs[-2]}.filtered_variants.csv'])
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
    dist = pdist(df, dist_nan)
    Z = linkage(np.nan_to_num(dist, 10), 'ward')
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


def get_ground_truth(in_files, c_max, min_cells):
    data = {}
    for in_file in in_files:
        if int(re.search('c(\d+)\.r\d+', in_file).group(1)) != c_max:
            continue
        if re.search('COMPASS', in_file):
            alg = 'COMPASS'
            try:
                dbt = re.search('c\d+\.r\d+\.d(0.\d+)', in_file).group(1)
            except AttributeError:
                continue
            assign, gt = get_COMPASS_results(in_file)
        elif re.search('BnpC', in_file):
            alg = 'BnpC'
            dbt = '-1'
            assign, gt = get_BnpC_results(in_file)
        elif re.search('SCG', in_file):
            alg = 'SCG'
            try:
                dbt = re.search('c\d+\.r\d+\.d(0.\d+)', in_file).group(1)
            except AttributeError:
                continue
            assign, gt = get_SCG_results(in_file)

        filter_small_clusters(assign, gt, min_cells)

        if not alg in data:
            data[alg] = {}
        if not dbt in data[alg]:
            data[alg][dbt] = {'assign': [], 'geno': []}

        data[alg][dbt]['assign'].append(assign)
        data[alg][dbt]['geno'].append(gt)

    gr_truth = {}
    for alg, dbt_data in data.items():
        gr_truth[alg] = {}
        for dbt, rel_data in dbt_data.items():
            gr_truth[alg][dbt] = {
                'cluster': get_coclustering_matrix(rel_data['assign']),
                'geno': get_geno_mean(rel_data['geno']),
                'assign': rel_data['assign']
            }

    return gr_truth


def filter_small_clusters(assign, gt, min_cells):
    cl_all, cl_size = np.unique(assign.values.flatten(), return_counts=True)
    for cl_i, size_i in zip(cl_all, cl_size):
        if size_i < min_cells:
            if assign.ndim == 2:
                rm_cells = (assign[0] == cl_i) & (assign[1] == 0)
                sgt_cells = ~(rm_cells) & (np.sum(assign == cl_i, axis=1) == 1)

                assign[rm_cells] = np.nan
                assign.replace(cl_i, 0, inplace=True)
                gt[rm_cells] = np.nan
                for cell, cell_data in assign[sgt_cells].iterrows():
                    if cell_data[0] == 0:
                        assign.loc[cell] = [cell_data[1], 0]
                    sgt_cl = assign.loc[cell, 0]
                    gt.loc[cell] = (gt[assign[0] == sgt_cl].iloc[0]).values
            else:
                rm_cells = (assign == cl_i)
                assign[rm_cells] = np.nan
                gt[rm_cells] = np.nan


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
            rm_cells = np.argwhere(np.sum(~np.isfinite(assign), axis=1)).flatten()
            vals.append((assign[None, :,0] == assign[:, None, 0]).astype(float) \
                + (assign[None, :,0] == assign[:, None, 1]).astype(float))
        else:
            rm_cells = np.argwhere(~np.isfinite(assign)).flatten()
            vals.append((assign[None,:] == assign[:,None]).astype(float))
        if rm_cells.size > 1:
            vals[-1][rm_cells] = np.nan
            vals[-1][:,rm_cells] = np.nan


    return pd.DataFrame(np.nanmean(np.stack(vals), axis=0), index=cells, columns=cells)


def dist_nan(u, v, scale=True):
    valid = ~np.isnan(u) & ~np.isnan(v)
    if valid.sum() == 0:
        return np.nan
    if scale:
        return euclidean(u[valid], v[valid]) / np.sqrt(np.sum(2**2 * valid.size))
    else:
        return euclidean(u[valid], v[valid])


def match_gt(true_df, pred_df):
    cells = pred_df.index
    SNPs = set(pred_df.columns) & set(true_df.columns)
    dist = dist_nan(
        true_df.loc[cells, SNPs].values.ravel(), pred_df.loc[:,SNPs].values.ravel())
    return 1 - dist, len(SNPs)


def match_assign(true_df, pred_df):
    cells = pred_df.index
    pwcm = get_coclustering_matrix([pred_df])
    dist = dist_nan(true_df.loc[cells, cells].values.ravel(), pwcm.values.ravel())
    return 1 - dist


def match_metric(metric, true_df, pred_df):
    scores = []
    for gt_assign in true_df:
        val_true = gt_assign.dropna().index
        val_red = pred_df.dropna().index
        rel = set(val_true) & set(val_red)

        try:
            _ = pred_df.shape[1]
        except IndexError:
            x = pred_df
            y = gt_assign
        else:
            # import pdb; pdb.set_trace()
            # pred_df.loc[pred_df[1] == 0, 1] = pred_df.loc[pred_df[1] == 0, 0]
            # gt_assign.loc[gt_assign[1] == 0, 1] = gt_assign.loc[gt_assign[1] == 0, 0]

            x = pred_df.where(lambda x: x[1] == 0, 0)[0]
            y = gt_assign.where(lambda x: x[1] == 0, 0)[0]

        try:
            scores.append(metric(
                y.loc[rel].values.ravel(), x.loc[rel].values.ravel()))
        except:
            import pdb; pdb.set_trace()

    return np.mean(scores)


def main(in_files, out_file, min_cells_all, plot_pwcm=False):
    c_max = max([int(re.search('/c(\d+)\.r\d+', i).group(1)) for i in in_files])

    gr_tr_all = {}
    data = []
    for min_cells in min_cells_all:
        gr_tr = get_ground_truth(in_files, c_max, min_cells)
        gr_tr_all[min_cells] = gr_tr

        for i, in_file in tqdm(enumerate(sorted(in_files))):
            if re.search('COMPASS', in_file):
                alg = 'COMPASS'
                assign, gt = get_COMPASS_results(in_file)
                try:
                    cell_no, run_no, dbt = re.search(
                        '/c(\d+)\.r(\d+)\.d(0.\d+)_cell', in_file).groups()
                except AttributeError:
                    continue
            elif re.search('BnpC', in_file):
                alg = 'BnpC'
                assign, gt = get_BnpC_results(in_file)
                cell_no, run_no = re.search('BnpC/c(\d+)\.r(\d+)/', in_file).groups()
                dbt = '-1'
            elif re.search('SCG', in_file):
                alg = 'SCG'
                assign, gt = get_SCG_results(in_file)
                try:
                    cell_no, run_no, dbt = re.search(
                        'SCG/c(\d+)\.r(\d+)\.d(0.\d+)/', in_file).groups()
                except AttributeError:
                    continue

            filter_small_clusters(assign, gt, min_cells)

            cl_no = np.isfinite(np.unique(assign)).sum()
            if alg == 'COMPASS':
                cl_no -= 1

            gt_score, SNP_overlap = match_gt(gr_tr[alg][dbt]['geno'], gt)
            assign_score = match_assign(gr_tr[alg][dbt]['cluster'], assign)
            AR_score = match_metric(adjusted_rand_score,
                gr_tr[alg][dbt]['assign'], assign)
            Mi_score = match_metric(adjusted_mutual_info_score,
                gr_tr[alg][dbt]['assign'], assign)
            data.append([int(cell_no), int(run_no), min_cells, alg, float(dbt),
                cl_no, gt_score, SNP_overlap, assign_score, AR_score, Mi_score])

    df = pd.DataFrame(data,
        columns=['cells', 'run', 'min. #cells per cluster', 'algorithm',
            'doublet rate', '#Cluster', 'Genotyping Sim.', '#SNPs overlapping',
            'Clustering Sim.', 'ARI', 'Mutual info'])
    df.sort_values(['cells', 'run', 'min. #cells per cluster', 'algorithm'],
        inplace=True)

    if not out_file:
        out_dir_raw = next(i for i in in_files if '/BnpC/' in i)
        out_dir = out_dir_raw.split('/BnpC/')[0]
        out_file = os.path.join(out_dir, 'summary.tsv')
    df.to_csv(out_file, sep='\t', index=False)

    if plot_pwcm:
        out_dir = os.path.dirname(out_file)
        for min_cells, gr_tr in gr_tr_all.items():
            for alg, dbt_data in gr_tr.items():
                for dbt, rel_data in dbt_data.items():
                    out_file_cocl = os.path.join(out_dir,
                        f'grTr_coclustering_minClSize{min_cells}_{alg}_d{dbt}.png')
                    plot_coclustering_matrix(rel_data['cluster'], out_file_cocl)


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
    parser.add_argument('-o', '--out_file', type=str, default='',
        help='Output file. Default = <INPUT>/summary.tsv')
    parser.add_argument('-m', '--min_cells', type=int, default=[1], nargs='+',
        help='Minimum number of cells to consider a cluster. Default = [1].')
    parser.add_argument('-p', '--plot_groundTruth', action='store_true',
        help='Plot pairwise coclustering matrices.')
    return parser.parse_args()


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.log:
            import sys
            sys.stderr = open(snakemake.log[0], 'w')
        main(snakemake.input, snakemake.output[0], snakemake.params.min_cells)
    else:
        args = parse_args()
        main(get_files(args.input), args.out_file, args.min_cells,
            args.plot_groundTruth)