#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy.spatial.distance import cityblock
from sklearn.metrics.cluster import adjusted_rand_score, adjusted_mutual_info_score


BNPC_PATTERN = r'BnpC/c(\d+)\.r(\d+)(?:\.relevant)*/'
COMPASS_PATTERN = r'COMPASS/c(\d+)\.r(\d+)(?:\.relevant)*\.d(0.\d+)_cell'
SCG_PATTERN = r'SCG/c(\d+)\.r(\d+)(?:\.relevant)*\.d(0.\d+)/'


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
    # Read cell assignments
    cells_att = pd.read_csv(assign_file, sep='\t', index_col=0)
    cells_att.sort_index(inplace=True)
    cells_att['node'] = cells_att \
        .apply(lambda x: float(x.node) if x.doublet == 'no' else np.nan, axis=1)

    # Read node genotypes
    gt_file = assign_file.replace('cellAssignments', 'nodes_genotypes')
    gt_nodes = pd.read_csv(gt_file, sep='\t')
    gt_nodes.set_index(gt_nodes['node'].str.replace('Node ', '').astype(int),
        inplace=True)
    gt_nodes.drop('node', axis=1, inplace=True)
    gt_nodes.loc[np.nan] = np.nan

    gt = gt_nodes.loc[cells_att['node']].set_index(cells_att.index)
    
    # Rename COMPASS output to match SNP annotation
    var_file = re.sub(r'\.d(\d\.)+.*$', '.filtered_variants.csv',
        assign_file.replace('COMPASS', 'samples'))
    df_var = pd.read_csv(var_file, dtype='str')
    snps = df_var['CHR'] + ':' + df_var['POS'] + ':' + df_var['REF'] + '/' + df_var['ALT']
    compass_annot = df_var.apply(lambda x: f'chr{x["CHR"]}_{x['POS']}({x["REGION"]})', axis=1)
    
    gt.rename(dict(zip(compass_annot, snps)), axis=1, inplace=True)
    gt = gt[sorted(gt.columns)]

    return cells_att['node'].sort_index(), gt.sort_index()


def get_BnpC_results(assign_file):
    chain = 'mean'
    estimator = 'posterior'
    # Read cell assignment
    df = pd.read_csv(assign_file, sep='\t', index_col=[0,1])
    assign = np.array(df.loc[(chain, estimator), 'Assignment'].split(' '),
        dtype=int)

    # Read cell genotype
    gt_file = os.path.join(os.path.dirname(assign_file),
        f'genotypes_{estimator}_{chain}.tsv')
    gt = pd.read_csv(gt_file, sep='\t', index_col=0).T
    gt = gt[sorted(gt.columns)]

    # Read SNPs
    path_dirs = os.path.abspath(assign_file).split(os.sep)
    sample_file = os.path.join(*[os.sep] + path_dirs[:-3] \
        + ['samples', f'{path_dirs[-2]}.filtered_variants.csv'])
    # Read cell names
    with open(sample_file, 'r') as f:
        cells = f.readline().strip().split(',')[7:]
    gt.index = cells

    return pd.Series(assign, index=cells).sort_index(), gt.sort_index()


def get_SCG_results(assign_file):
    # Read cell assignment
    cells_att = pd.read_csv(assign_file, sep='\t', index_col=0)
    cells_att.sort_index(inplace=True)
    cells_att.rename({i: i.strip('\'b') for i in cells_att.index}, axis=0, inplace=True)

    
    # Read cell genotype
    gt_file = os.path.join(os.path.dirname(assign_file), 'Cluster_genotypes.tsv')
    gt = pd.read_csv(gt_file, sep='\t', index_col=0)
    gt.rename({i: i.strip('\'b') for i in gt.columns}, axis=1, inplace=True)
    gt = gt[sorted(gt.columns)]
    gt = gt.loc[cells_att.loc[:, 'cluster_id']]

    return cells_att['cluster_id'].sort_index(), gt.set_index(cells_att.index).sort_index()


def get_ground_truth(in_files, c_max):
    data = {}
    for in_file in in_files:
        # Keep only max cell numbers
        if int(re.search(r'c(\d+)\.r\d+', in_file).group(1)) != c_max:
            continue

        if re.search('COMPASS', in_file):
            alg = 'COMPASS'
            try:
                dbt = re.search(COMPASS_PATTERN, in_file).group(3)
            except AttributeError:
                continue
            run_no = re.search(COMPASS_PATTERN, in_file).group(2)
            assign, gt = get_COMPASS_results(in_file)
        elif re.search('BnpC', in_file):
            alg = 'BnpC'
            dbt = '-1'
            run_no = re.search(BNPC_PATTERN, in_file).group(2)
            assign, gt = get_BnpC_results(in_file)
        elif re.search('SCG', in_file):
            alg = 'SCG'
            try:
                dbt = re.search(SCG_PATTERN, in_file).group(3)
            except AttributeError:
                continue
            run_no = re.search(SCG_PATTERN, in_file).group(2)
            assign, gt = get_SCG_results(in_file)

        assign.name = f'{alg}_{run_no}'
        gt.name = f'{alg}_{run_no}'
        
        if '.relevant' in in_file:
            subset = 'relevant'
            SNPs_relevant = set(gt.columns)
        else:
            subset = 'all'
            SNPs_all = set(gt.columns)
        
        if not alg in data:
            data[alg] = {}
        if not subset in data[alg]:
            data[alg][subset] = {}
        if not dbt in data[alg][subset]:
            data[alg][subset][dbt] = {'assign': [], 'gt': []}
        data[alg][subset][dbt]['assign'].append(assign)
        data[alg][subset][dbt]['gt'].append(gt)

    # Take cells from last read gt (sorted + same in all files)
    cells = gt.index.values
    # Merge data
    for alg, alg_data in data.items():
        for subset, subset_data in alg_data.items():
            for dbt, dbt_data in subset_data.items():
                data[alg][subset][dbt]['assign'] = pd.concat(dbt_data['assign'], axis=1)
                data[alg][subset][dbt]['gt'] = np.stack([i.values for i in dbt_data['gt']])
        
    return data, (SNPs_all, SNPs_relevant), cells



def _filter_small_cl(df, min_cells):
    cl_ids, cl_sizes = np.unique(df.values, return_counts=True)
    val_cells = ~df.isna()
    for cl_i, size_i in zip(cl_ids, cl_sizes):
        if size_i < min_cells:
            cl_cells = df[df == cl_i].index.values
            val_cells.loc[cl_cells] = False
    return val_cells


def match_metric(metric, true_dfs, pred_df, min_cells):
    val_cells_pred = _filter_small_cl(pred_df, min_cells)
    scores = []
    for run_no, true_df in true_dfs.T.iterrows():
        val_cells_true = _filter_small_cl(true_df, min_cells)
        val_cells = val_cells_true & val_cells_pred
        scores.append(metric(true_df.loc[val_cells], pred_df.loc[val_cells]))

    return np.mean(scores).round(4)


def get_hamming_dist(true_arrs, pred_df, SNPs_all, cells_all):
    # Get index (ground truth data) of overlapping cells and SNPs with subsample
    # Dropping cells called as doublets in COMPASS
    cells_idx = np.where(np.isin(cells_all, pred_df.dropna().index))[0]
    SNPs_idx = np.where(np.isin(list(SNPs_all), pred_df.columns))[0]
    # Drop SNPs that were not called in the ground truth data
    SNPS_pred_idx = np.isin(pred_df.columns, list(SNPs_all))
    pred_arr = pred_df.loc[:, SNPS_pred_idx].dropna().values
    
    scores = []
    for true_arr in true_arrs[:,cells_idx,:][:,:,SNPs_idx]:
        dbt_true = np.isnan(true_arr)
        # For COMPASS: Exclude cells that are called as doublets
        if dbt_true.any():
            dbt_true_cells = np.argwhere(dbt_true.sum(axis=1) == SNPs_idx.size).ravel()
            true_vals = np.delete(true_arr, dbt_true_cells, axis=0).ravel()
            pred_vals = np.delete(pred_arr, dbt_true_cells, axis=0).ravel()
        else:
            true_vals = true_arr.ravel()
            pred_vals = pred_arr.ravel()
        
        scores.append(cityblock(true_vals, pred_vals) / pred_vals.size)
    
    return np.mean(scores).round(4)


def main(in_files, out_file, min_cells_all, plot_pwcm=False):
    c_max = max([int(re.search(r'/c(\d+)\.r\d+', i).group(1)) for i in in_files])
    data_gt, SNPs_gt_subsets, cells_all = get_ground_truth(in_files, c_max)

    data = []
    for i, in_file in enumerate(sorted(in_files)):
        if re.search('COMPASS', in_file):
            alg = 'COMPASS'
            assign, gt = get_COMPASS_results(in_file)
            cell_no, run_no, dbt = re.search(COMPASS_PATTERN, in_file).groups()
        elif re.search('BnpC', in_file):
            alg = 'BnpC'
            assign, gt = get_BnpC_results(in_file)
            cell_no, run_no = re.search(BNPC_PATTERN, in_file).groups()
            dbt = '-1'
        elif re.search('SCG', in_file):
            alg = 'SCG'
            assign, gt = get_SCG_results(in_file)
            cell_no, run_no, dbt = re.search(SCG_PATTERN, in_file).groups()

        if '.relevant' in in_file:
            subset = 'relevant'
            SNPs_gt = SNPs_gt_subsets[1]
        else:
            subset = 'all'
            SNPs_gt = SNPs_gt_subsets[0]
        
        assign.name = f'{alg}_{run_no}'
        gt.name = f'{alg}_{run_no}'

        SNPs = set(gt.columns)
        SNP_overlap = len(SNPs_gt & SNPs)
        SNP_missed = len(SNPs_gt - SNPs)
        SNP_new = len(SNPs - SNPs_gt)
        run_data = [int(cell_no), subset, int(run_no), alg, float(dbt),
            SNP_overlap, SNP_missed, SNP_new]

        for min_cells in min_cells_all:
            cl_no = (assign.fillna(-1).value_counts().drop(-1, errors='ignore') \
                >= min_cells).sum()

            AR_score = match_metric(adjusted_rand_score,
                data_gt[alg][subset][dbt]['assign'], assign, min_cells)
            Mi_score = match_metric(adjusted_mutual_info_score,
                data_gt[alg][subset][dbt]['assign'], assign, min_cells)
            hamming_dist = get_hamming_dist(data_gt[alg][subset][dbt]['gt'], gt,
                SNPs_gt, cells_all)
            res_data = [min_cells, cl_no, AR_score, Mi_score, hamming_dist]
            data.append(run_data + res_data)
    run_cols = ['cells', 'subset', 'run', 'algorithm', 'doublet rate',
        'SNPs_overlap', 'SNPs_missed', 'SNPs_new']
    res_cols = ['min. #cells per cluster', 'Cluster', 'ARI', 'Mutual info',
        'Norm. Hamming distance']
    df = pd.DataFrame(data, columns=run_cols + res_cols)
    df.sort_values(['cells', 'subset', 'run', 'min. #cells per cluster',
        'algorithm'], inplace=True)

    if not out_file:
        out_dir_raw = next(i for i in in_files if '/BnpC/' in i)
        out_dir = out_dir_raw.split('/BnpC/')[0]
        out_file = os.path.join(out_dir, 'summary.tsv')
    df.to_csv(out_file, sep='\t', index=False)


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
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input directory (output dir of sampling run) or all run files.')
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
        if len(args.input) == 1 and os.path.isdir(args.input[0]):
            in_files = get_files(args.input[0])
        else:
            in_files = args.input
        main(in_files, args.out_file, args.min_cells, args.plot_groundTruth)
