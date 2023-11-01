#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set_context('talk')
from matplotlib.colors import LinearSegmentedColormap


DISTS = ['Mismatch', 'Manhatten', 'Euclidean']

SGT = ['0->0', '1->1', '2->2']
SGT_WRONG = ['0->1', '0->2', '1->0', '1->2', '2->0', '2->1']
SGT_DBT = ['0->0+1', '0->1+2', '0->0+2', '1->0+1', '1->1+2', '1->0+2',
    '2->0+1', '2->1+2', '2->0+2']

DBT = ['0+1->0+1', '0+2->0+2', '1+2->1+2']
DBT_WRONG = ['0+1->0+2', '0+1->1+2', '0+2->0+1', '0+2->1+2', '1+2->0+1', '1+2->0+2']
DBT_SGT = ['0+1->0', '0+1->1', '0+1->2', '0+2->0', '0+2->1', '0+2->2',
    '1+2->0', '1+2->1', '1+2->2']


def get_scSplit_results(df):
    pred = df['Cluster'].apply(lambda x: x.split('-')[-1]).values
    dbt = [i.split('-')[1] for i in df['Cluster'].unique() if 'DBL' in i][0]
    true_cl = np.array(
        ['+'.join(re.findall('.pat(\d+)', i)) for i in df.index.values])
    true_dbt = np.char.find(true_cl, '+') == 1

    pat_map = {}
    for cl_id in np.unique(pred):
        true_pats = true_cl[np.argwhere(pred == cl_id).flatten()]
        pats, cnts = np.unique(true_pats, return_counts=True)

        sugg_pat = pats[np.argmax(cnts)]
        if sugg_pat not in pat_map.values() or pats.size == 1:
            pat_map[cl_id] = sugg_pat
        else:
            sugg_frac = np.max(cnts) / cnts.sum()
            old_cl = [i for i,j  in pat_map.items() if j == sugg_pat][0]
            old_pats, old_cnts = np.unique(
                true_cl[np.argwhere(pred == old_cl).flatten()],
                return_counts=True)
            old_frac = np.max(old_cnts) / old_cnts.sum()
            if old_frac >= sugg_frac:
                pat_map[cl_id] = pats[np.argsort(cnts)[-2]]
            else:
                pat_map[old_cl] = old_pats[np.argsort(old_cnts)[-2]]
                pat_map[cl_id] = sugg_pat
        if '+' in pat_map[cl_id]:
            pat_map[cl_id] = 'DBT'

    right = 0
    wrong = 0
    wrong_sgt = 0
    wrong_dbt = 0
    res = {}
    for i, j in pat_map.items():
        cells = np.argwhere(pred == i).flatten()
        if i == dbt:
            right += true_dbt[cells].sum()
            wrong += (~true_dbt[cells]).sum()
            wrong_sgt += (~true_dbt[cells]).sum()
        else:
            right += ((pred[cells] == i) & (true_cl[cells] == j)).sum()
            wrong += ((pred[cells] == i) & (true_cl[cells] != j)).sum()
            wrong_dbt += true_dbt[cells].sum()

        for cl, cnt in zip(*np.unique(true_cl[cells], return_counts=True)):
            res[f'{cl}->{pat_map[i]}'] = cnt

    return (right, wrong, wrong_sgt, wrong_dbt), res


def get_soup_results(df):
    pred = df['assignment'] \
        .apply(lambda x: x if not '/' in x else '+'.join(sorted(x.split('/')))) \
        .values
    dbt = (df['status'] == 'doublet').values
    true_cl = np.array(
        ['+'.join(sorted(re.findall('.pat(\d+)', i))) for i in df.index.values])
    true_dbt = np.char.find(true_cl, '+') == 1

    pat_map = {}
    for cl_id in np.unique(pred):
        true_pats = true_cl[np.argwhere(pred == cl_id).flatten()]
        pats, cnts = np.unique(true_pats, return_counts=True)

        sugg_pat = pats[np.argmax(cnts)]
        if sugg_pat not in pat_map.values() or pats.size == 1:
            pat_map[cl_id] = sugg_pat
        else:
            sugg_frac = np.max(cnts) / cnts.sum()
            old_cl = [i for i,j  in pat_map.items() if j == sugg_pat][0]
            old_pats, old_cnts = np.unique(
                true_cl[np.argwhere(pred == old_cl).flatten()],
                return_counts=True)
            old_frac = np.max(old_cnts) / old_cnts.sum()
            if old_frac >= sugg_frac:
                pat_map[cl_id] = pats[np.argsort(cnts)[-2]]
            else:
                pat_map[old_cl] = old_pats[np.argsort(old_cnts)[-2]]
                pat_map[cl_id] = sugg_pat

    right = 0
    wrong = 0
    wrong_sgt = 0
    wrong_dbt = 0
    res = {}
    for i, j in pat_map.items():
        cells = np.argwhere(pred == i).flatten()

        if '+' in j:
            right += true_dbt[cells].sum()
            wrong += (~true_dbt[cells]).sum()
            wrong_sgt += (~true_dbt[cells]).sum()
        else:
            right += ((pred[cells] == i) & (true_cl[cells] == j)).sum()
            wrong += ((pred[cells] == i) & (true_cl[cells] != j)).sum()
            wrong_dbt += true_dbt[cells].sum()

        for cl, cnt in zip(*np.unique(true_cl[cells], return_counts=True)):
            res[f'{cl}->{pat_map[i]}'] = cnt

    return (right, wrong, wrong_sgt, wrong_dbt), res



def get_vireo_results(df):
    df['best_singlet'] = df['best_singlet'].str.replace('donor', '')
    df['best_doublet'] = df['best_doublet'] \
        .apply(lambda x: '+'.join(sorted([i[-1] for i in x.split(',')])))

    sgt = ((df['donor_id'] != 'doublet') & (df['donor_id'] != 'unassigned')).values
    dbt = (df['donor_id'] == 'doublet').values

    pred = df['best_singlet'].where(sgt, df['best_doublet']).values
    # Add 'unassigned' predctions
    pred[(df['donor_id'] == 'unassigned').values] = '-1'
    true_cl = np.array(
        ['+'.join(sorted(re.findall('.pat(\d+)', i))) for i in df.index.values])
    true_dbt = np.char.find(true_cl, '+') == 1

    pat_map = {'-1': '-1'}
    for cl_id in np.unique(pred):
        if cl_id == '-1':
            continue
        true_pats = true_cl[np.argwhere(pred == cl_id).flatten()]
        pats, cnts = np.unique(true_pats, return_counts=True)

        sugg_pat = pats[np.argmax(cnts)]
        if sugg_pat not in pat_map.values() or pats.size == 1:
            pat_map[cl_id] = sugg_pat
        else:
            sugg_frac = np.max(cnts) / cnts.sum()
            old_cl = [i for i,j  in pat_map.items() if j == sugg_pat][0]
            old_pats, old_cnts = np.unique(
                true_cl[np.argwhere(pred == old_cl).flatten()],
                return_counts=True)
            old_frac = np.max(old_cnts) / old_cnts.sum()
            if old_frac >= sugg_frac:
                pat_map[cl_id] = pats[np.argsort(cnts)[-2]]
            else:
                pat_map[old_cl] = old_pats[np.argsort(old_cnts)[-2]]
                pat_map[cl_id] = sugg_pat

    right = 0
    wrong = 0
    wrong_sgt = 0
    wrong_dbt = 0
    res = {}

    for i, j in pat_map.items():
        cells = np.argwhere(pred == i).flatten()

        if '+' in j:
            right += true_dbt[cells].sum()
            wrong += (~true_dbt[cells]).sum()
            wrong_sgt += (~true_dbt[cells]).sum()
        elif j == '-1':
            wrong += cells.size
            wrong_sgt += (~true_dbt[cells]).sum()
            wrong_dbt += (true_dbt[cells]).sum()
        else:
            right += ((pred[cells] == i) & (true_cl[cells] == j)).sum()
            wrong += ((pred[cells] == i) & (true_cl[cells] != j)).sum()
            wrong_dbt += true_dbt[cells].sum()

        for cl, cnt in zip(*np.unique(true_cl[cells], return_counts=True)):
            res[f'{cl}->{pat_map[i]}'] = cnt

    return (right, wrong, wrong_sgt, wrong_dbt), res


def main(args):
    for in_file in sorted(args.input):
        df_new = pd.read_csv(in_file, sep='\t', index_col=0)
        df_new.rename({'Euclidean': 'demoTape'}, axis=1, inplace=True)
        df_new.drop('file', axis=1, inplace=True, errors='ignore')

        if 'scSplit' in in_file:
            rep = int(re.search('/rep(\d+)/', in_file).group(1))
            scSplit_short, scSplit_long  = get_scSplit_results(df_new)
            df_sum_new = pd.DataFrame(scSplit_short).T.rename({0: 'scSplit'})
            df_new = pd.DataFrame([scSplit_long]).T.rename({0: 'scSplit'}, axis=1)
            df_new['rep'] = rep
            df_new['algorithm'] = 'scSplit'
        elif 'souporcell' in in_file:
            rep = int(re.search('/rep(\d+)/', in_file).group(1))
            soup_short, soup_long  = get_soup_results(df_new)
            df_sum_new = pd.DataFrame(soup_short).T.rename({0: 'souporcell'})
            df_new = pd.DataFrame([soup_long]).T.rename({0: 'souporcell'}, axis=1)
            df_new['rep'] = rep
            df_new['algorithm'] = 'souporcell'
        elif 'vireo' in in_file:
            rep = int(re.search('/rep(\d+)/', in_file).group(1))
            vireo_short, vireo_long  = get_vireo_results(df_new)
            df_sum_new = pd.DataFrame(vireo_short).T.rename({0: 'vireo'})
            df_new = pd.DataFrame([vireo_long]).T.rename({0: 'vireo'}, axis=1)
            df_new['rep'] = rep
            df_new['algorithm'] = 'vireo'
        else:
            rep = int(re.search('rep(\d+).distance', in_file).group(1))
            right = df_new.loc[SGT].sum() + df_new.loc[DBT].sum() + df_new.loc[DBT_WRONG].sum()
            wrong = df_new.loc[SGT_WRONG].sum() + df_new.loc[SGT_DBT].sum() \
                + df_new.loc[DBT_SGT].sum()
            wrong_sgt = df_new.loc[SGT_DBT].sum()
            wrong_dbt = df_new.loc[DBT_SGT].sum()
            df_sum_new = pd.concat([right, wrong, wrong_sgt, wrong_dbt], axis=1)
            df_new['rep'] = rep
            df_new['algorithm'] = 'demoTape'

        df_sum_new['rep'] = rep

        df_sum_new = df_sum_new.reset_index().set_index(['rep', 'index'])
        df_sum_new.columns = ['right', 'wrong', 'sgt2dbt', 'dbt2sgt']

        try:
            df = df.append(df_new, sort=False)
            df_sum = df_sum.append(df_sum_new, sort=False)
        except NameError:
            df_sum = df_sum_new
            df = df_new

    if not args.output:
        args.output = os.path.join(os.path.dirname(args.input[0]), 'summary.tsv')

    df_sum.to_csv(args.output, sep='\t', index=True)
    out_full_base, out_full_end = os.path.splitext(args.output)
    df.fillna(0).to_csv(f'{out_full_base}_full{out_full_end}', sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input files from demultplexing runs.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for similarity df. default = <INPUT[0]_DIR>.summary.tsv')
    parser.add_argument('-op', '--output_plot', action='store_true',
        help='Output file for heatmap with dendrogram to "<INPUT>.hm.png".')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)