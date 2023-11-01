#!/usr/bin/env python3

import argparse
import copy
from itertools import combinations
import os
import re

import numpy as np
import matplotlib
# matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import pdist, euclidean
from scipy.special import comb
import seaborn as sns
from seaborn import clustermap


COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#E4761A' , '#2FBF85'] # '#A4CA55'
COLORS = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
    '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']
COLORS_STR = ['red', 'blue', 'green', 'orange', 'mint']

WHITELIST = []
WL_IDS = []


FILE_EXT = 'png'
FONTSIZE = 30
DPI = 300
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'lines.linewidth': 2,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'large',
        'xtick.labelsize': 'medium',
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 50
})
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.bottom'] = True


def dist_nan(u, v, scale=True):
    valid = ~np.isnan(u) & ~np.isnan(v)
    if scale:
        return euclidean(u[valid], v[valid]) / np.sqrt(np.sum(2**2 * valid.size))
    else:
        return euclidean(u[valid], v[valid])


def merge_gt(gt_in):
    # 00=0, 11=1, 22=2,
    # 33=3; 03=0, 13=1, 23=2;
    # 01=1, 02=1, 12=1
    if gt_in[0] == gt_in[1]:
        return gt_in[0]
    elif np.isnan(gt_in[0]):
        return gt_in[1]
    elif np.isnan(gt_in[1]):
        return gt_in[0]
    else:
        return 1


def main(args):
    in_files = []

    if len(args.input) == 1 and os.path.isdir(args.input[0]):
        for file in os.listdir(args.input[0]):
            if file.endswith('_variants.csv'):
                in_files.append(os.path.join(args.input[0], file))
    else:
        in_files = args.input

    for file_name in in_files:
        df_new = dist_clustering(file_name, args)
        df_new['file'] = file_name
        try:
            df = pd.concat([df, df_new])
        except NameError:
            df = df_new

    if args.output_summary:
        df.to_csv(args.output_summary, sep='\t')


def dist_clustering(in_file, args):
    df = pd.read_csv(in_file, index_col=[0, 1], dtype={'CHR': str})

    labels = df.apply(lambda x:
            #f'{x.name[0]}_{x.name[1]} ({x["NAME"]}): {x["REF"]}>{x["ALT"]}',
            f'chr{x.name[0]}:{x.name[1]} {x["REF"]}>{x["ALT"]}',
        axis=1)
    if WHITELIST:
        var_ids = df.apply(lambda x: f'{x.name[0]}_{x.name[1]}', axis=1)
        WL_IDS.extend(np.argwhere(np.isin(var_ids.values, WHITELIST)).flatten())

    df.drop(['REF', 'ALT', 'REGION', 'NAME', 'FREQ'], axis=1, inplace=True)

    split_pat = '[\.,+]'
    true_cl = []
    for cell in df.columns:
        if 'pat' in cell:
            s1 = int(re.split('[\.,]', cell.split('+')[0])[-1][3:])
        else:
            s1 = 0

        if '+' in cell:
            s2 = int(re.split('[\.,]', cell.split('+')[1])[-1][3:])
            true_cl.append('+'.join([str(j) for j in sorted([s1, s2])]))
        else:
            true_cl.append(s1)
    true_cl = np.array(true_cl)

    gt = df.applymap(lambda x: float(x.split(':')[-1])).values.T
    gt[gt == 3] = np.nan

    dist = pdist(gt, dist_nan)
    Z = linkage(dist, 'ward')

    # i.e., if synthetic data without doublets
    if args.no_doublets:
        cl_no = args.clusters
    else:
        dbt_cl_no = int(comb(args.clusters, 2))
        cl_no = args.clusters + dbt_cl_no

    while True:
        done = True
        clusters = cut_tree(Z, n_clusters=cl_no).flatten()
        # plot_heatmap(gt, Z, true_cl, [i for i in map(str, clusters)], labels, None)

        profiles = np.zeros(shape=(cl_no, gt.shape[1]))
        for cl_id in range(cl_no):
            profiles[cl_id] = np.nanmean(gt[clusters == cl_id], axis=0)

        dbt_profiles = np.zeros(shape=(int(comb(cl_no, 2)), gt.shape[1]))
        dbt_combs = []

        for i, combo in enumerate(combinations(range(cl_no), 2)):
            dbt_profiles[i] = np.apply_along_axis(merge_gt, axis=0,
                arr=profiles[combo,:].round())
            dbt_combs.append(f'{combo[0]}+{combo[1]}')

        dbt_cl_map = {}
        dbt_ids = []
        if not args.no_doublets and args.clusters > 1:
            df_dbt = pd.DataFrame([], columns=dbt_combs, index=range(cl_no))
            for i in range(cl_no):
                df_dbt.loc[i] = np.apply_along_axis(dist_nan, 1, dbt_profiles,
                    profiles[i])

            for i in range(cl_no):
                for j in dbt_combs:
                    if str(i) in j:
                        df_dbt.loc[i, j] = np.nan

            for i in dbt_combs:
                i1, i2 = map(int, i.split('+'))
                if dist_nan(profiles[i1], profiles[i2]) < 0.1:
                    df_dbt.loc[:, i] = np.nan

            # Run until (3, 2) doublet identified or no more value to consider
            while len(dbt_ids) < comb(args.clusters, 2):
                # Check if dist array is empty all all na
                if df_dbt.size == 0 or df_dbt.isna().sum().sum() == df_dbt.size:
                    print('Cannot solve doublet clusters')
                    done = False
                    break

                cl_idx, dbt_idx = np.where(df_dbt == df_dbt.min().min())
                cl_id = df_dbt.index[cl_idx[0]]
                dbt_id = df_dbt.columns[dbt_idx[0]]

                dcl1, dcl2 = map(int, dbt_id.split('+'))

                hom_dcl1 = set(np.argwhere(profiles[dcl1] > 1.5).flatten())
                hom_dcl2 = set(np.argwhere(profiles[dcl2] > 1.5).flatten())

                hom12_id = np.array(list(hom_dcl1 & hom_dcl2))
                hom1_id = np.array(list(hom_dcl1 - hom_dcl2))
                hom2_id = np.array(list(hom_dcl2 - hom_dcl1))
                # Check if both hom match
                if hom12_id.size > 0:
                    hom12_match = profiles[cl_id][hom12_id] > 1.5
                else:
                    hom12_match = []
                # Check if hom on cl1 match
                if hom1_id.size > 0:
                    hom1_match = (profiles[cl_id][hom1_id] >= 0.5) \
                        & (profiles[cl_id][hom1_id] < 1.5)
                else:
                    hom1_match = []
                # Check if hom on cl2 match
                if hom2_id.size > 0:
                    hom2_match = (profiles[cl_id][hom2_id] >= 0.5) \
                        & (profiles[cl_id][hom2_id] < 1.5)
                else:
                    hom2_match = []

                hom_match = np.concatenate([hom12_match, hom1_match, hom2_match]).mean()
                print(f'Doublet cluster similarity: {hom_match:.3f}')
                if hom_match < 0.95:
                    print(f'Hom. match between {cl_id} and {dbt_id} low: {hom_match:.3f}')
                    done = False
                    break

                dbt_ids.append(cl_id)
                dbt_cl_map[cl_id] = dbt_id
                # Remove cluster row
                rm_rows = [cl_id, dcl1, dcl2]
                df_dbt.drop(rm_rows, axis=0, inplace=True, errors='ignore')
                # Remove doublet rows including cluster
                rm_cols = [dbt_id] + [i for i in dbt_combs if str(cl_id) in i]
                df_dbt.drop(rm_cols, axis=1, inplace=True, errors='ignore') # Ignore error on already dropped labels

        if done:
            break
        cl_no += 1
        print(f'Increasing cuttree clusters to {cl_no}')

    def get_sgt_dist_matrix(ids, profiles):
        mat = np.zeros((ids.size, ids.size))
        for i, j in enumerate(ids):
            mat[i] = np.apply_along_axis(dist_nan, 1, profiles[ids], profiles[j])
        mat[np.tril_indices(ids.size)] = np.nan
        return mat

    def update_dbt_map(dbt_cl_map, new_cl, del_cl):
        dbt_cl_map_new = {}
        consistent = True
        for new_key, val in dbt_cl_map.items():
            if new_key > del_cl:
                new_key -= 1

            sgt1, sgt2 = map(int, val.split('+'))
            if sgt1 == del_cl:
                sgt1 = new_cl
            elif sgt1 > del_cl:
                sgt1 = sgt1 - 1
            if sgt2 == del_cl:
                sgt2 = new_cl
            elif sgt2 > del_cl:
                sgt2 = sgt2 - 1

            # By merging, a doublet cluster is not a doublet anymore
            if sgt1 == sgt2:
                consistent = False
            else:
                dbt_cl_map_new[new_key] = '+'.join(sorted([str(sgt1), str(sgt2)]))

        return dbt_cl_map_new, consistent


    dbt_ids = np.array(dbt_ids)
    sgt_ids = np.array([i for i in np.arange(cl_no) if not i in dbt_ids])
    # Merge surplus single clusters
    while sgt_ids.size > args.clusters:
        sgt_dists = get_sgt_dist_matrix(sgt_ids, profiles)
        merge_failed = True
        while True:
            if np.nansum(sgt_dists) == 0:
                merge_failed = True
                sgt_dists = get_sgt_dist_matrix(sgt_ids, profiles)

            sgt1_idx, sgt2_idx = np.where(sgt_dists == np.nanmin(sgt_dists))
            cl1 = sgt_ids[sgt1_idx[0]]
            cl2 = sgt_ids[sgt2_idx[0]]

            dbt_cl_map_new, consistent = update_dbt_map(dbt_cl_map, cl1, cl2)
            if len(dbt_cl_map_new.values()) != len(set(dbt_cl_map_new.values())):
                consistent = False

            if not consistent and not merge_failed:
                sgt_dists[sgt1_idx, sgt2_idx] = np.nan
                continue
            break

        dbt_cl_map = dbt_cl_map_new

        # Remove cl2
        clusters[clusters == cl2] = cl1
        clusters[clusters > cl2] -= 1

        sgt_ids = np.delete(sgt_ids, sgt2_idx[0])
        profiles = np.delete(profiles, cl2, axis=0)

        dbt_ids[dbt_ids > cl2] -= 1
        sgt_ids[sgt_ids > cl2] -= 1

        if dbt_ids.size != len(dbt_cl_map):
            print('Removing 1 obsolete doublet cluster')
            rem_dbt_cl = (set(dbt_ids) - set(dbt_cl_map.keys())).pop()

            clusters[clusters == rem_dbt_cl] = cl1
            clusters[clusters > rem_dbt_cl] -= 1
            dbt_ids = np.delete(dbt_ids, np.argwhere(dbt_ids == rem_dbt_cl))
            sgt_ids[sgt_ids > rem_dbt_cl] -= 1
            dbt_ids[dbt_ids > rem_dbt_cl] -= 1

            dbt_cl_map, _ = update_dbt_map(dbt_cl_map, cl1, rem_dbt_cl)

        if len(dbt_cl_map.values()) != len(set(dbt_cl_map.values())):
            print('Removing 1 equal doublet cluster')
            eq_dbt = [i for i,j in dbt_cl_map.items() \
                if list(dbt_cl_map.values()).count(j) > 1]

            clusters[clusters == eq_dbt[1]] = eq_dbt[0]
            clusters[clusters > eq_dbt[1]] -= 1
            dbt_ids = np.delete(dbt_ids, np.argwhere(dbt_ids == eq_dbt[1]))
            sgt_ids[sgt_ids > eq_dbt[1]] -= 1
            dbt_ids[dbt_ids > eq_dbt[1]] -= 1

            dbt_cl_map.pop(eq_dbt[1])
            dbt_cl_map, _ = update_dbt_map(dbt_cl_map, eq_dbt[0], eq_dbt[1])

            profiles[eq_dbt[0]] = np.nanmean(gt[clusters == eq_dbt[0]], axis=0)

        # Update profile
        profiles[cl1] = np.nanmean(gt[clusters == cl1], axis=0)

    # Rename clusters to range from 0 to n
    cl_map = {j: str(i) for i, j in enumerate(sgt_ids)}
    for i, j in dbt_cl_map.items():
        cl1, cl2 = map(int, j.split('+'))
        cl_map[i] = '+'.join(sorted([cl_map[int(cl1)], cl_map[int(cl2)]]))

    clusters_str = clusters.astype('<U3')
    for i, j in enumerate(clusters):
        clusters_str[i] = cl_map[j]

    # Print cluster summaries to stdout
    for cl_id, cl_size in zip(*np.unique(clusters, return_counts=True)):
        cl_name = cl_map[cl_id]
        if '+' in cl_name:
            cl1, cl2 = cl_name.split('+')
            cl_color = f'{COLORS_STR[int(cl1)]}+{COLORS_STR[int(cl2)]}'
        else:
            cl_color = COLORS_STR[int(cl_name)]
        print(f'Cluster {cl_name} ({cl_color}): {cl_size: >4} cells ' \
            f'({cl_size / df.shape[1] * 100: >2.0f}%)')

        gt_cl_called = np.isfinite(gt[clusters == cl_id]).sum(axis=0)
        for geno in [0, 1, 2]:
            gt_cl = gt[clusters == cl_id] == geno
            print(f'\tGT {geno} - Avg./cell {gt_cl.sum(axis=1).mean(): >4.1f}, '
                f'# 95% clonal: {(gt_cl.sum(axis=0) > gt_cl_called * 0.95).sum()}')

    pat_map = {}
    # Snythetic multiplexed data: ground truth known
    if np.unique(true_cl).size > 1:
        for cl_id in np.unique(clusters):
            true_pats = true_cl[np.argwhere(clusters == cl_id).flatten()]
            pats, cnts = np.unique(true_pats, return_counts=True)

            sugg_pat = pats[np.argmax(cnts)]
            if sugg_pat not in pat_map.values() or pats.size == 1:
                pat_map[cl_id] = sugg_pat
            else:
                sugg_frac = np.max(cnts) / cnts.sum()
                old_cl = [i for i,j  in pat_map.items() if j == sugg_pat][0]
                old_pats, old_cnts = np.unique(
                    true_cl[np.argwhere(clusters == old_cl).flatten()],
                    return_counts=True)
                old_frac = np.max(old_cnts) / old_cnts.sum()
                if old_frac >= sugg_frac:
                    pat_map[cl_id] = pats[np.argsort(cnts)[-2]]
                else:
                    pat_map[old_cl] = old_pats[np.argsort(old_cnts)[-2]]
                    pat_map[cl_id] = sugg_pat

    if np.unique(true_cl).size > 1: # works only for snythetic multiplexed data
        res_cell = {}
        if args.clusters == 2:
            pat_strs = ['0', '1', '0+1']
        elif args.clusters == 3:
            pat_strs = ['0', '1', '2', '0+1', '0+2', '1+2']
        elif args.clusters == 4:
            pat_strs = ['0', '1', '2', '3',
                '0+1', '0+2', '0+3', '1+2', '1+3', '2+3']
        else:
            raise RuntimeError('Not implemented for >3 clusters')
        counts = {f'{i}->{j}': 0 for i in pat_strs for j in pat_strs}

        for cell_id, cl_id in enumerate(clusters):
            cell_true = true_cl[cell_id]

            if cl_id in dbt_ids:
                cl_dist = np.apply_along_axis(dist_nan, 1, profiles[sgt_ids],
                    profiles[cl_id])
                dbt_cl = sgt_ids[np.argsort(cl_dist)[:2]]
                pred_cl = sorted([pat_map[dbt_cl[0]], pat_map[dbt_cl[1]]])
                pred = f'{pred_cl[0]}+{pred_cl[1]}'
            else:
                pred = pat_map[cl_id]
            try:
                counts[f'{cell_true}->{pred}'] += 1
            except KeyError:
                counts[f'{cell_true}->{pred}'] = 1
        res_cell['Euclidean'] = counts

        res_cell_df = pd.DataFrame(res_cell)
        res_cell_df.index.name = 'Assignment'
        res_cell_df.columns.name = 'Distance'
    else:
        res_cell_df = pd.DataFrame([])

    if cl_no > 1:
        # Print assignment to tsv file
        cell_str = '\t'.join([i for i in df.columns.values])
        cl_str = '\t'.join([i for i in clusters_str])
        with open(f'{in_file}.assignments.tsv', 'w') as f:
            f.write(f'Barcode\t{cell_str}\nCluster\t{cl_str}')

        # Safe SNV profiles to identy patients
        idx = [f'{cl_map[i]} ({COLORS_STR[int(cl_map[i])]})' for i in sgt_ids]
        pd.DataFrame(profiles[sgt_ids], index=idx, columns=labels).round(2) \
            .to_csv(f'{in_file}.profiles.tsv', sep='\t')

    if args.output_plot or args.show_plot:
        if args.output_plot:
            out_file = f'{in_file}.heatmap.{FILE_EXT}'
        else:
            out_file = None
        snv_order, cell_order = plot_heatmap(gt, Z, true_cl, clusters_str,
            labels, out_file)

        avg_idx = []
        for cl_id in np.unique(clusters):
            avg_idx.append(
                (cl_id,
                np.argwhere(clusters_str[cell_order] == cl_map[cl_id]).mean())
            )
        row_order = [i[0] for i in sorted(avg_idx, key=lambda x: x[1])]
        if WL_IDS:
            profile_gt = profiles.round()[row_order][:,WL_IDS][:,snv_order]
        else:
            profile_gt = profiles.round()[row_order][:,snv_order]
        plot_heatmap(profile_gt, None, np.zeros(cl_no),
            [cl_map[i] for i in row_order], labels[snv_order],
            f'{in_file}.heatmap.profiles.{FILE_EXT}', cluster=False)

    return res_cell_df


def get_row_cols(assignment):
    row_colors1 = []
    row_colors2 = []
    for i in assignment:
        if '+' in i:
            s1, s2 = [int(j) for j in i.split('+')]
            if s1 > s2:
                row_colors1.append(COLORS[s2])
                row_colors2.append(COLORS[s1])
            else:
                row_colors1.append(COLORS[s1])
                row_colors2.append(COLORS[s2])
        else:
            row_colors1.append(COLORS[int(i)])
            row_colors2.append('#FFFFFF')
    return row_colors1, row_colors2


def plot_heatmap(gt, Z, true_cl, clusters_str, labels, out_file, cluster=True):
    myColors = ('#EAEAEA', '#fed976', '#fc4e2a', '#800026')
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    cl_no = np.unique(clusters_str).size
    # Get row color colors based on assignment
    if cl_no > 1:
        r_colors = []
        if np.unique(true_cl).size > 1:
            r_colors.extend(get_row_cols(true_cl))
        r_colors.extend(get_row_cols(clusters_str))
    else:
        r_colors = None

    if WL_IDS and gt.shape[1] != len(WL_IDS):
        df_plot = np.nan_to_num(gt[:, WL_IDS], nan=-1)
        labels_plot = labels[WL_IDS]
    else:
        df_plot = np.nan_to_num(gt, nan=-1)
        labels_plot = labels

    cm = clustermap(
        df_plot,
        row_linkage=Z,
        row_cluster=cluster,
        col_cluster=cluster,
        vmin=-1, vmax=2,
        row_colors=r_colors,
        cmap=cmap,
        figsize=(25, 10),
        xticklabels=labels_plot,
        cbar_kws={'shrink': 0.5, 'drawedges': True}
    )
    cm.ax_heatmap.set_ylabel('Cells')
    if df_plot.shape[0] > 50:
        cm.ax_heatmap.set_yticks([])
    cm.ax_heatmap.set_xlabel('SNVs')

    cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xticklabels(),
        rotation=45, fontsize=5, ha='right', va='top')

    cm.ax_col_dendrogram.set_visible(False)

    colorbar = cm.ax_heatmap.collections[0].colorbar
    colorbar.set_ticks([-0.65, 0.15, 0.85, 1.65])
    colorbar.set_ticklabels([r' $-$', '0|0', '0|1', '1|1'])

    cm.fig.tight_layout()
    if out_file:
        print(f'Saving heatmap to: {out_file}')
        cm.fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()

    if cluster:
        return cm.dendrogram_col.reordered_ind, cm.dendrogram_row.reordered_ind


def update_whitelist(wl_file):
    df = pd.read_csv(wl_file, sep='\t')

    if df.size == 0:
        if df.shape[1] > 1:
            for snv in df.columns:
                WHITELIST.append('_'.join(snv.split('_')[:2]))
        else:
            for snv in df.columns.values[0].split(';'):
                WHITELIST.append('_'.join(snv.split('_')[:2]))
    elif df.shape[1] == 1:
        df = pd.read_csv(wl_file, sep=',')
        ids = df['CHR'].astype(str) + '_' + df['POS'].astype(str)
        WHITELIST.extend(ids)
    else:
        for _, (sample, snvs) in df.iterrows():
            if '+' in sample:
                continue
            for snv in snvs.split(';'):
                WHITELIST.append('_'.join(snv.split('_')[:2]))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input _variants.csv file(s).')
    parser.add_argument('-n', '--clusters', type=int, default=1,
        help='Number of clusters to define. Default = 1.')
    parser.add_argument('-wl', '--whitelist', type=str, default='',
        help='Whitelist containint SNVs for plotting. Default = None.')
    parser.add_argument('-nd', '--no_doublets', action='store_true',
        help='If set, no doublets are assumed.')

    output = parser.add_argument_group('output')
    output.add_argument('-op', '--output_plot', action='store_true',
        help='Output file for heatmap with dendrogram to "<INPUT>.hm.png".')
    output.add_argument('-sp', '--show_plot', action='store_true',
        help='Show heatmap with dendrogram at stdout.')
    output.add_argument('-os', '--output_summary', type=str,
        help='Output file for demultiplexing summary (simualted data only).')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.whitelist:
        update_whitelist(args.whitelist)
    main(args)