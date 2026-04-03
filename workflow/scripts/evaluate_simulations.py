#!/usr/bin/env python3

import argparse
import os
import re

import numpy as np
import pandas as pd
from sklearn.metrics.cluster import v_measure_score



def replace_doublets(arr, value='-1'):
    return np.where(np.char.find(arr.astype(str), '+') != -1,
        value, arr.astype(str))


def obtain_assignment_results(pred_all, true_all):
    # Set doublets to -1 cluster
    true_all = replace_doublets(true_all)
    # DoubletD special case: doublet cl = -1, singlet = 0
    if 'singlet' in pred_all:
        true_all = np.where(true_all == '-1', '-1', '0')
        pred_all = np.where(pred_all == 'doublet', '-1', '0')    
    else:
        pred_all = replace_doublets(pred_all)
    
    return v_measure_score(true_all, pred_all)


def get_predictions_vireo(df):
    df['best_singlet'] = df['best_singlet'].str.replace('donor', '')
    df['best_doublet'] = df['best_doublet'] \
        .apply(lambda x: '+'.join(sorted([i[-1] for i in x.split(',')])))

    sgt = ((df['donor_id'] != 'doublet') & (df['donor_id'] != 'unassigned')).values
    dbt = (df['donor_id'] == 'doublet').values

    pred_all = df['best_singlet'].where(sgt, df['best_doublet']).values
    # Add 'unassigned' predictions
    pred_all[(df['donor_id'] == 'unassigned').values] = '-2'
    true_all = np.array(
        ['+'.join(sorted(re.findall(r'.pat(\d+)', i))) for i in df.index.values])
    return pred_all, true_all


def get_predictions_demoTape(df):
    pred_all = df.loc['Cluster'].values.astype(str)
    true_all = np.array(
        ['+'.join(sorted(re.findall(r'.pat(\d+)', i))) for i in df.columns.values])
    return pred_all, true_all


def get_predictions_scSplit(df):
    # Wrong axis pd.concat
    if sum([i.startswith('Cluster') for i in df.columns]) > 1:
        data = []
        for i in range(0, df.shape[1] + 1, 2):
            data.append(df.reset_index().iloc[:,i:i+2].dropna()
                .rename(columns=lambda x: x.split('.')[0]))
        df = pd.concat(data).set_index('Barcode')
    pred_all = df['Cluster'].apply(lambda x: x.split('-')[-1]).values
    true_all = np.array(
        ['+'.join(re.findall(r'.pat(\d+)', i)) for i in df.index.values])
    return pred_all, true_all


def get_predictions_soup(df):
    pred_all = df['assignment'] \
        .apply(lambda x: x if not '/' in x else '+'.join(sorted(x.split('/')))) \
        .values
    true_all = np.array(
        ['+'.join(sorted(re.findall(r'.pat(\d+)', i))) for i in df.index.values])
    return pred_all, true_all

def get_predictions_doubletD(df):
    pred_all = df['prediction'].values
    true_all = np.array(
        ['+'.join(sorted(re.findall(r'.pat(\d+)', i))) for i in df.index.values])
    return pred_all, true_all


def main(in_files, out_file):
    for in_file in sorted(in_files):
        df_new = pd.read_csv(in_file, sep='\t', index_col=0)
        df_new.rename({'Euclidean': 'demoTape'}, axis=1, inplace=True)
        df_new.drop('file', axis=1, inplace=True, errors='ignore')

        if 'scSplit' in in_file:
            name = 'scSplit'
            pred_all, true_all = get_predictions_scSplit(df_new)
        elif 'souporcell' in in_file:
            name = 'souporcell'
            pred_all, true_all = get_predictions_soup(df_new)
        elif 'vireo' in in_file:
            name = 'vireo'
            pred_all, true_all = get_predictions_vireo(df_new)
        elif 'doubletD' in in_file:
            name = 'doubletD'
            pred_all, true_all = get_predictions_doubletD(df_new)
        else:
            name = 'demoTape'
            pred_all, true_all = get_predictions_demoTape(df_new)
        try:
            rep = int(re.search(r'/rep(\d+)', in_file).group(1))
        except:
            rep = -1
        res = obtain_assignment_results(pred_all, true_all)

        res_new = pd.DataFrame({'rep': rep, 'algorithm': name, 'v_measure': res},
            index=[0])

        try:
            df = pd.concat([df, res_new], axis=0, ignore_index=True)
        except NameError:
            df = res_new

    if not out_file:
        out_file = os.path.join(os.path.dirname(in_files[0]), 'summary.tsv')

    df.to_csv(out_file, sep='\t', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input files from demultplexing runs.'),
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file for results df. default = <INPUT[0]_DIR>.summary.tsv')

    return parser.parse_args()


if __name__ == '__main__':
    if 'snakemake' in globals():
        if snakemake.log:
            import sys
            sys.stderr = open(snakemake.log[0], 'w')
        main(snakemake.input, snakemake.output[0])
    else:
        args = parse_args()
        main(args.input, args.output)