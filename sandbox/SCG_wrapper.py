#!/usr/bin/env python3

import argparse
from datetime import datetime
import os
import subprocess
import sys

import yaml
import h5py
import numpy as np
import pandas as pd

# ------------------------------------------------------------------------------
# INIT AND OUTPUT FUNCTIONS
# ------------------------------------------------------------------------------

def preprocess_SCG(args):
    if not args.output:
        res_dir = f'{datetime.now():%Y%m%d_%H:%M:%S}_SCG'
        args.output = os.path.join(os.path.dirname(args.input), res_dir)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    data_out_file = os.path.join(args.output, f'{os.path.basename(args.input)}.gz')
    data = pd.read_csv(args.input, sep=',', index_col=0)

    data.T.to_csv(data_out_file, index_label='cell_id', sep='\t')

    cfg = {
        'num_clusters': args.clusters,
        'alpha_prior': [(1 - args.doublets) / args.doublets, 1], # Beta dist: E[X] = a / (a + b); b = 1
        'kappa_prior': 1,
        'data': {
            'snv': {
                'file': data_out_file,
                'gamma_prior':
                    [[98, 1, 1],
                    [25, 50, 25],
                    [1, 1, 98]]
                ,
                'state_prior': [1, 1, 1]
            }
        }
    }

    if args.doublets:
        cfg['model'] = 'doublet'
        cfg['state_map'] = {'snv': {
            0: [[0, 0]],
            1: [[0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1]],
            2: [[2, 2]]}
        }

    cfg_file = os.path.join(args.output, 'config.yaml')
    with open(cfg_file, 'w') as f:
        yaml.dump(cfg, f)
    args.config_file = cfg_file


def run_SCG(args):
    hdf_file = os.path.join(args.output, 'scg_run.hdf')
    cmmd = f'scg fit --in-file {args.config_file} --out-file {hdf_file} ' \
        f'--max-iters {args.steps}'
    if not args.silent:
        print(f'output directory:\n{args.output}')
        print(f'\nShell command:\n{cmmd}\n')

    scg = subprocess.Popen(cmmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = scg.communicate()
    scg.wait()

    if not stdout or str(stdout).split('\\n')[-2] != 'Converged':
        print('\nSCG didnt converge:')
        [print(i) for i in str(stderr).split('\\n')]
        raise RuntimeError('SCG Error')

    if not args.silent:
        print('\nSCG stdout:')
        for i in str(stdout).split('\\n'):
            print(i)


def postprocess_SCG(out_dir):
    hdf_file = os.path.join(out_dir, 'scg_run.hdf')
    cell_df, cluster_df = load_results_dfs(hdf_file)

    # safe assignment
    assign_file = os.path.join(out_dir, 'assignments.tsv')
    cell_df.to_csv(assign_file, sep='\t', index=False, float_format='%.4f')

    # Load genotypes
    df_geno = pd.DataFrame()
    for cluster_no, cluster_data in cluster_df.groupby('cluster_id'):
        cluster_genotype = cluster_data \
            .loc[cluster_data.groupby('event_id')['genotype_prob'].idxmax()] \
            .set_index('event_id')
        df_geno[cluster_no] = cluster_genotype['genotype']

    df_geno.T.to_csv(os.path.join(out_dir, 'Cluster_genotypes.tsv'), '\t',
        index_label='Cluster')


def main(args):
    if args.results:
        postprocess_SCG(args.results)
    else:
        preprocess_SCG(args)
        run_SCG(args)
        postprocess_SCG(args.output)


def load_results_dfs(file_name):
    cell_df = load_cell_df(file_name)
    cluster_df = load_cluster_df(file_name)

    cluster_map = dict(zip(
        cell_df['cluster_id'].unique(), np.arange(cell_df['cluster_id'].nunique())
    ))

    cell_df["cluster_id"] = cell_df["cluster_id"].map(cluster_map)
    cell_df = cell_df.sort_values(by="cluster_id")

    cluster_df = cluster_df[cluster_df["cluster_id"].isin(cluster_map.keys())]
    cluster_df["cluster_id"] = cluster_df["cluster_id"].map(cluster_map)
    cluster_df = cluster_df.sort_values(by=["cluster_id", "event_id"])

    return cell_df, cluster_df


def load_cell_df(file_name):
    with h5py.File(file_name, 'r') as fh:
        model = fh["meta"].attrs["model"]
        cell_ids = fh["/data/cell_ids"][()]

        Z = fh["/var_params/Z"][()]
        if model == "doublet":
            K = fh["meta"].attrs["K"]
            Y = fh["/var_params/Y"][()]
            # Keep only the columns for single clusters
            Z = Z[:,:K]
            Z = Z / Z.sum(axis=1)[:, np.newaxis]

            df = pd.DataFrame({
                "cluster_id": Z.argmax(axis=1),
                "cluster_prob": Z.max(axis=1),
                "doublet_prob": Y[1]
            })
        else:
            df = pd.DataFrame({
                "cluster_id": Z.argmax(axis=1),
                "cluster_prob": Z.max(axis=1)
            })

        df.insert(0, "cell_id", cell_ids)

    return df


def load_cluster_df(file_name):
    with h5py.File(file_name, 'r') as fh:
        data_types = list(fh["/data/event_ids"].keys())

        df = []
        for dt in data_types:
            event_ids = fh[f"/data/event_ids/{dt}"][()]

            G = fh[f"/var_params/G/{dt}"][()]
            # MAP assignments
            G_idx = G.argmax(axis=0)
            G_idx = pd.DataFrame(G_idx, columns=event_ids)
            G_idx = G_idx.stack().reset_index()
            G_idx.columns = "cluster_id", "event_id", "genotype"

            # Probs
            G_prob = G.max(axis=0)
            G_prob = pd.DataFrame(G_prob, columns=event_ids)
            G_prob = G_prob.stack().reset_index()
            G_prob.columns = "cluster_id", "event_id", "genotype_prob"

            # Merge
            dt_df = pd.merge(G_idx, G_prob, on=["cluster_id", "event_id"])
            dt_df["data_type"] = dt
            df.append(dt_df)

    return pd.concat(df)


def parse_args():
    parser = argparse.ArgumentParser(
          description='*** Wrapper for SCG (single-cell genotyper). ***'
    )
    parser.add_argument('input', help='Absolute or relative path to input data. ' \
           'Input data is a n x m matrix (n = cells, m = mutations) with 1|0, ' \
           'representing whether a mutation is present in a cell or not. Matrix ' \
           'elements need to be separated by a whitespace or tabulator.'
    )
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory. Default = "<DATA_DIR>/<TIMESTAMP>_SCG".'
    )
    parser.add_argument('-k', '--clusters', type=int, default=50,
        help='Maximum number of clusters. Default = 50.'
    )
    parser.add_argument('-db', '--doublets', type=float, default=0.08,
        help='Expected doublet rate. Default = 0.08'
    )
    parser.add_argument('-s', '--steps', type=int, default=10000,
        help='Maximum number of ELBO optimization iterations. Default = 10000.'
    )
    parser.add_argument('-si', '--silent', action='store_true',
        help='Print status massages to stdout. Default = True'
    )
    parser.add_argument('-r', '--results', type=str, default='',
        help='Only postprocess results dir.'
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    main(args)