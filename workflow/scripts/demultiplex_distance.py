#!/usr/bin/env python3

import argparse
from itertools import combinations
import os
import re

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from scipy.cluster.hierarchy import linkage, cut_tree, leaves_list
from scipy.spatial.distance import pdist, euclidean, cityblock
from scipy.special import comb
from scipy.stats import binom, kstest
import seaborn as sns
from seaborn import clustermap

EPSILON = np.finfo(np.float64).resolution

COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#E4761A' , '#2FBF85'] # '#A4CA55'
COLORS = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
    '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']

COLORS = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
COLORS_STR = ['red', 'blue', 'green', 'orange', 'mint']


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


class demoTape:
    def __init__(self, in_file, cl_no, distance='manhatten'):
        self.sgt_cl_no = cl_no
        self.dbt_cl_no = int(comb(self.sgt_cl_no, 2))

        # Get full data
        df = pd.read_csv(in_file, index_col=[0, 1], dtype={'CHR': str})
        self.SNPs = df.apply(
            lambda x: f'chr{x.name[0]}:{x.name[1]} {x["REF"]}>{x["ALT"]}',
            axis=1)
        self.rel_SNPs = np.ones(self.SNPs.size, dtype=bool)
        df.drop(['REF', 'ALT', 'REGION', 'NAME', 'FREQ'], axis=1, inplace=True)
        self.cells = df.columns.values

        # Get true clusters, if simulated data
        true_cl = []
        for cell in df.columns:
            if 'pat' in cell:
                s1 = int(re.split(r'[\.,]', cell.split('+')[0])[-1][3:])
            else:
                s1 = 0

            if '+' in cell:
                s2 = int(re.split(r'[\.,]', cell.split('+')[1])[-1][3:])
                true_cl.append('+'.join([str(j) for j in sorted([s1, s2])]))
            else:
                true_cl.append(s1)
        self.true_cl = np.array(true_cl)

        # Init variables
        self.Z = np.array([])
        self.assignment = np.array([])
        self.profiles = np.array([])
        self.sgt_ids = np.array([])
        self.dbt_ids = np.array([])
        self.dbt_map = {}

        self.get_data(df)


    def __str__(self):
        out_str = '\ndemoTape:\n' \
            f'\tFile: {self.in_file}\n' \
            f'\t# Samples: {self.cl_no}\n' \
            f'\tCells: {self.cells}:\n' \
            f'\tSNPs: {self.SNPs}\n'
        return out_str


    @staticmethod
    def metric(c1, c2):
        pass


    def get_data(self, df):
        pass


    def get_pairwise_dists(self):
        pass


    def get_profile(self, cl):
        pass


    def check_hom_match(self, scl, dcl):
        pass



    def demultiplex(self):
        self.init_dendrogram()
        self.identify_doublets()
        self.merge_surplus_singlets()


# -------------------------------- DENDROGRAM ----------------------------------

    def init_dendrogram(self):
        dist = self.get_pairwise_dists()
        self.Z = np.clip(linkage(np.nan_to_num(dist, 0.5), method='ward'), 0, None)


# ---------------------------- IDENTIFY DOUBLETS -------------------------------

    def identify_doublets(self):
        cl_no = self.sgt_cl_no + self.dbt_cl_no
        # Generate clusters until (n, 2) doublets clusters are identified
        while True:
            self.set_assignment(cl_no)
            self.set_dbt_ids()

            if self.dbt_ids.size == self.dbt_cl_no:
                break
            elif cl_no == (self.sgt_cl_no + self.dbt_cl_no) * 3:
                print(f'Could not identify all doublets.')
                self.set_assignment(self.sgt_cl_no + self.dbt_cl_no)
                self.set_dbt_ids()
                break
            cl_no += 1
            print(f'Increasing cuttree clusters to {cl_no}')


    def set_assignment(self, cl_no):
        self.assignment = cut_tree(self.Z, n_clusters=cl_no).flatten()
        clusters = np.unique(self.assignment)
        self.profiles = np.zeros(shape=(clusters.size, self.prof_len))
        for cl_id, cl in enumerate(clusters):
            self.profiles[cl_id] = self.get_profile(cl)



    def set_dbt_ids(self):
        cl_size  = np.unique(self.assignment, return_counts=True)[1] / self.cells.size
        cl_no = cl_size.size

        dbt_map = {}
        dbt_ids = []

        dbt_profiles = np.zeros(shape=(int(comb(cl_no, 2)), self.prof_len))
        dbt_combs = []

        for i, combo in enumerate(combinations(range(cl_no), 2)):
            dbt_profiles[i] = self.get_profile(combo)
            dbt_combs.append(combo)

        df_dbt = pd.DataFrame([], columns=dbt_combs, index=range(cl_no))
        for i in range(cl_no):
            df_dbt.loc[i] = np.apply_along_axis(self.metric, 1, dbt_profiles,
                self.profiles[i])
            for j in dbt_combs:
                # set doublet combo dist to np.nan if:
                # 1. Singlet is included in Dbt cluster
                # 2. Dbt cluster size is larger than one of the two clusters
                if i in j or  cl_size[i] > cl_size[[*j]].min() * 2:
                    df_dbt.loc[i, [j]] = np.nan

        # Set dbt combo dists to np.nan if their profiles are very close
        for i in dbt_combs:
            if self.metric(self.profiles[i[0]], self.profiles[i[1]]) < 0.05:
                df_dbt.loc[:, [i]] = np.nan

        # Run until (n, 2) doublet identified or no more value to consider
        while len(dbt_ids) < self.dbt_cl_no:
            # Check if dist array is empty or all np.nan
            if df_dbt.size == 0 or df_dbt.isna().sum().sum() == df_dbt.size:
                print('Cannot solve doublet clusters')
                break

            cl_idx, dbt_idx = np.where(df_dbt == df_dbt.min().min())
            cl_id = df_dbt.index[cl_idx[0]]
            dbt_id = df_dbt.columns[dbt_idx[0]]

            # Allow max 1 missmatch
            if self.check_hom_match(cl_id, dbt_id) > 1:
                df_dbt.loc[cl_id, [dbt_id]] = np.nan
                continue

            dbt_ids.append(cl_id)
            dbt_map[cl_id] = dbt_id
            # Remove cluster row
            rm_rows = [cl_id] + list(dbt_id)
            df_dbt.drop(rm_rows, axis=0, inplace=True, errors='ignore')
            # Remove doublet rows including cluster
            rm_cols = [dbt_id] + [i for i in dbt_combs if cl_id in i]
            df_dbt.drop(rm_cols, axis=1, inplace=True, errors='ignore') # Ignore error on already dropped labels

        self.sgt_ids = np.array([i for i in np.unique(self.assignment) \
            if not i in dbt_ids])
        self.dbt_ids = np.array(dbt_ids)
        self.dbt_map = dbt_map

        self.plot_heatmap(rel_SNPs_only=True)
        import pdb; pdb.set_trace()


    def get_cl_map(self):
        if self.sgt_ids.size == 0:
            cl_map = {i: i for i in np.unique(self.assignment)}
            assignment_str = [str(i) for i in self.assignment]
        else:
            # Rename clusters to range from 0 to n
            cl_map = {j: str(i) for i, j in enumerate(self.sgt_ids)}
            for i, j in self.dbt_map.items():
                cl_map[i] = '+'.join(sorted([cl_map[j[0]], cl_map[j[1]]]))

            assignment_str = np.full(self.cells.size, '-1', dtype='<U3')
            for i, j in enumerate(self.assignment):
                assignment_str[i] = cl_map[j]
        return cl_map, assignment_str



# ------------------------- MERGE SURPLUS SINGLETS -----------------------------

    def merge_surplus_singlets(self):
        # Merge surplus single clusters
        while self.sgt_ids.size > self.sgt_cl_no:
            sgt_dists = self.get_sgt_dist_matrix()
            while True:
                # Connect solve doublet merging without conflict
                if np.nansum(sgt_dists) == 0:
                    sgt_dists = self.get_sgt_dist_matrix()
                    sgt1_idx, sgt2_idx = np.where(sgt_dists == np.nanmin(sgt_dists))
                    cl1 = self.sgt_ids[sgt1_idx[0]]
                    cl2 = self.sgt_ids[sgt2_idx[0]]

                    dbt_map_new, _ = self.update_dbt_map(cl1, cl2)
                # Possible clusters to merge left
                else:
                    sgt1_idx, sgt2_idx = np.where(sgt_dists == np.nanmin(sgt_dists))
                    cl1 = self.sgt_ids[sgt1_idx[0]]
                    cl2 = self.sgt_ids[sgt2_idx[0]]

                    dbt_map_new, consistent = self.update_dbt_map(cl1, cl2)
                    if len(dbt_map_new.values()) != len(set(dbt_map_new.values())):
                        consistent = False

                    if not consistent:
                        sgt_dists[sgt1_idx, sgt2_idx] = np.nan
                        continue

                break

            # if not consistent:
            #     self.plot_heatmap()
            #     import pdb; pdb.set_trace()

            self.dbt_map = dbt_map_new

            # Remove cl2
            self.assignment[self.assignment == cl2] = cl1
            self.assignment[self.assignment > cl2] -= 1

            self.sgt_ids = self.sgt_ids[np.where(self.sgt_ids != cl2)]
            self.profiles = np.delete(self.profiles, cl2, axis=0)

            self.dbt_ids[self.dbt_ids > cl2] -= 1
            self.sgt_ids[self.sgt_ids > cl2] -= 1

            if np.unique(self.dbt_ids).size != self.dbt_ids.size:
                print('Removing 1 obsolete doublet cluster')
                self.dbt_ids = np.unique(self.dbt_ids)
            elif self.dbt_ids.size != len(self.dbt_map):
                print('Removing 1 obsolete doublet cluster')
                rem_dbt_cl = (set(self.dbt_ids) - set(self.dbt_map.keys())).pop()

                self.assignment[self.assignment == rem_dbt_cl] = cl1
                self.assignment[self.assignment > rem_dbt_cl] -= 1
                self.dbt_ids = np.delete(self.dbt_ids,
                    np.argwhere(self.dbt_ids == rem_dbt_cl))
                self.sgt_ids[self.sgt_ids > rem_dbt_cl] -= 1
                self.dbt_ids[self.dbt_ids > rem_dbt_cl] -= 1

                self.dbt_map, _ = self.update_dbt_map(cl1, rem_dbt_cl)

            if len(self.dbt_map.values()) != len(set(self.dbt_map.values())):
                print('Removing 1 equal doublet cluster')
                print(self.dbt_map)
                eq_dbt = sorted([i for i,j in self.dbt_map.items() \
                    if list(self.dbt_map.values()).count(j) > 1])

                self.assignment[self.assignment == eq_dbt[1]] = eq_dbt[0]
                self.assignment[self.assignment > eq_dbt[1]] -= 1
                self.dbt_ids = np.delete(self.dbt_ids,
                    np.argwhere(self.dbt_ids == eq_dbt[1]))
                self.sgt_ids[self.sgt_ids > eq_dbt[1]] -= 1
                self.dbt_ids[self.dbt_ids > eq_dbt[1]] -= 1

                self.dbt_map.pop(eq_dbt[1])
                self.dbt_map, _ = self.update_dbt_map(eq_dbt[0], eq_dbt[1])

                cl1 = eq_dbt[0]

            # Update profile
            self.profiles[cl1] = self.get_profile(cl1)



    def get_sgt_dist_matrix(self):
        mat = np.zeros((self.sgt_ids.size, self.sgt_ids.size))
        for i, j in enumerate(self.sgt_ids):
            mat[i] = np.apply_along_axis(self.metric, 1,
                self.profiles[self.sgt_ids], self.profiles[j])
        mat[np.tril_indices(self.sgt_ids.size)] = np.nan
        return mat


    def update_dbt_map(self, new_cl, del_cl):
        dbt_map_new = {}
        consistent = True
        for new_key, val in self.dbt_map.items():
            if new_key > del_cl:
                new_key -= 1

            sgt1, sgt2 = val
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
                dbt_map_new[new_key] = tuple(sorted([sgt1, sgt2]))

        return dbt_map_new, consistent


    # ------------------------ HEATMAP GENERATION ------------------------------

    @staticmethod
    def get_cmap():
        pass


    def get_hm_data(self):
        pass


    @staticmethod
    def get_cm_specifics():
        pass


    @staticmethod
    def apply_cm_specifics(cm):
        pass


    def plot_heatmap(self, out_file='', cluster=True, rel_SNPs_only=False):
        cmap = self.get_cmap()

        # Get row color colors based on assignment
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

        r_colors = []
        if np.unique(self.true_cl).size > 1:
            r_colors.extend(get_row_cols(self.true_cl))

        _, assignment_str = self.get_cl_map()
        r_colors.extend(get_row_cols(assignment_str))

        df_plot, mask = self.get_hm_data()

        if self.Z.size > 0:
            Z = self.Z
        else:
            Z = None

        if rel_SNPs_only:
            df_plot = df_plot[:, self.rel_SNPs]
            mask = mask[:, self.rel_SNPs]
            x_ticks = self.SNPs[self.rel_SNPs]
        else:
            x_ticks = self.SNPs

        cm = clustermap(
            df_plot,
            row_linkage=self.Z,
            row_cluster=cluster,
            col_cluster=False,
            mask=mask,
            row_colors=r_colors,
            cmap=cmap,
            figsize=(25, 10),
            xticklabels=x_ticks,
            **self.get_cm_specifics()
        )

        cm.ax_row_colors.set_xlabel('Clusters', rotation=45, fontsize=10,
            ha='right', va='top')

        cm.ax_heatmap.set_facecolor('#5B566C')
        cm.ax_heatmap.set_ylabel('Cells')
        if df_plot.shape[0] > 50:
            cm.ax_heatmap.set_yticks([])
        cm.ax_heatmap.set_xlabel('SNPs')

        cm.ax_heatmap.set_xticklabels(cm.ax_heatmap.get_xticklabels(),
            rotation=45, fontsize=5, ha='right', va='top')

        cm.ax_col_dendrogram.set_visible(False)
        self.apply_cm_specifics(cm)

        cm.fig.tight_layout()
        if out_file:
            print(f'Saving heatmap to: {out_file}')
            cm.fig.savefig(out_file, dpi=DPI)
        else:
            plt.show()


    # -------------------------- GENERATE OUTPUT -------------------------------

    def safe_results(self, output):
        self.print_summary()
        self.safe_profiles(output)
        # Print assignment to tsv file
        _, assignment_str = self.get_cl_map()
        cell_str = '\t'.join([i for i in self.cells])
        cl_str = '\t'.join([i for i in assignment_str])
        order_str = '\t'.join([str(i) for i in leaves_list(self.Z)])
        with open(f'{output}.assignments.tsv', 'w') as f:
            f.write(f'Barcode\t{cell_str}\nCluster\t{cl_str}\nOrder\t{order_str}')


    def safe_profiles(self, output):
        pass



    def print_summary(self):
        pass


# ------------------------------------------------------------------------------


class demoTape_reads(demoTape):
    def __init__(self, in_file, cl_no):
        super().__init__(in_file, cl_no)
        self.prof_len = self.SNPs.size * 3
        self.rel_SNPs = self.get_relevant_SNPS()


    def __str__(self):
        return str(super()) + f'Distance metric: reads'


    def get_relevant_SNPS(self):
        # Identify SNPs that are symmetric/normal distributed: likely germline + ADO
        # Not informative for clustering
        VAF_z = np.nan_to_num(self.VAF - np.nanmean(self.VAF, axis=0), 0)
        q_vals_z = [kstest(VAF_z[:,i], -1 * VAF_z[:,i]).pvalue * self.SNPs.size \
            for i in range(self.SNPs.size)]
        rel_SNPs = np.array([i < 0.05 for i in q_vals_z])

        # Identify SNPs that are on the same read in most/all cells and remove 
        #   one of them from clustering
        for chrom, chrom_SNPs in self.SNPs[rel_SNPs].groupby('CHR'):
            if chrom_SNPs.size < 2:
                continue
            pos = chrom_SNPs.index.get_level_values('POS').values
            # SNPs at position id[x] and [x - 1] are on the same read
            same_read_SNP = np.argwhere((pos[1:] - pos[:-1]) < 275).ravel() + 1
            if same_read_SNP.size == 0:
                continue

            for SNP2_chrom_id in same_read_SNP:
                # Get row index in full data
                SNP1 = self.SNPs.index.get_loc(chrom_SNPs.index[SNP2_chrom_id - 1])
                SNP2 = self.SNPs.index.get_loc(chrom_SNPs.index[SNP2_chrom_id])
                
                # Calculate euclidean distance between VAF profiles
                valid = ~np.isnan(self.VAF[:,SNP1]) & ~np.isnan(self.VAF[:,SNP2])
                VAF_dist = euclidean(self.VAF[:,SNP1][valid], self.VAF[:,SNP2][valid]) \
                    / np.sqrt(np.sum(2**2 * valid.sum()))
                # If VAF profiles are very similar: remove second SNP
                if VAF_dist < 0.05:
                    rel_SNPs[SNP2] = False

        return rel_SNPs


    def get_data(self, df):
        self.ref = df.applymap(lambda x: int(x.split(':')[0])).values.T
        self.alt = df.applymap(lambda x: int(x.split(':')[1])).values.T
        self.dp = self.ref + self.alt
        VAF = (self.alt + EPSILON) / (self.dp + EPSILON)
        self.VAF = np.clip(np.where(self.dp > 0, VAF, np.nan), EPSILON,
            1 - EPSILON)
        self.RAF = 1 - self.VAF
        self.norm_const = np.insert(
            np.arange(1, self.dp.max() * 2 + 1) \
                * np.log(np.arange(1, self.dp.max() * 2 + 1)),
            0, np.nan)
        self.reads = np.hstack([self.ref, self.alt, self.dp])


    def metric(self, c1, c2):
        r1 = np.reshape(c1, (3, -1))[:, self.rel_SNPs]
        r2 = np.reshape(c2, (3, -1))[:, self.rel_SNPs]

        valid = (r1[2] > 0) & (r2[2] > 0)
        r1 = r1[:,valid]
        r2 = r2[:,valid]

        p1 = np.clip(r1[0] / r1[2], EPSILON, 1 - EPSILON)
        p2 = np.clip(r2[0] / r2[2], EPSILON, 1 - EPSILON)
        dp_total = r1[2] + r2[2]
        p12 = np.clip((r1[0] + r2[0]) / (dp_total), EPSILON, 1 - EPSILON)
        p12_inv = 1 - p12

        logl = r1[0] * np.log(p1 / p12) + r1[1] * np.log((1 - p1) / p12_inv) \
            + r2[0] * np.log(p2 / p12) + r2[1] * np.log((1 - p2) / p12_inv)

        norm = np.log(dp_total) * (dp_total) \
            - r1[2] * np.log(r1[2]) - r2[2] * np.log(r2[2])

        return np.sum(logl / norm) / valid.sum()


    def get_pairwise_dists(self):
        dist = []
        for i in np.arange(self.cells.size - 1):
            valid = (self.dp[i, self.rel_SNPs] > 0) & (self.dp[i+1:, self.rel_SNPs] > 0)
            dp_total = self.dp[i, self.rel_SNPs] + self.dp[i+1:, self.rel_SNPs]
            p12 = np.clip(
                (self.alt[i, self.rel_SNPs] + self.alt[i+1:, self.rel_SNPs] \
                    + EPSILON) / (dp_total + EPSILON),
                EPSILON, 1 - EPSILON)
            p12_inv = 1 - p12

            logl = self.alt[i, self.rel_SNPs] * np.log(self.VAF[i, self.rel_SNPs] / p12) \
                + self.ref[i, self.rel_SNPs] * np.log(self.RAF[i, self.rel_SNPs] / p12_inv) \
                + self.alt[i+1:, self.rel_SNPs] * np.log(self.VAF[i+1:, self.rel_SNPs] / p12) \
                + self.ref[i+1:, self.rel_SNPs] * np.log(self.RAF[i+1:, self.rel_SNPs] / p12_inv)

            norm = self.norm_const[dp_total] \
                - self.norm_const[self.dp[i, self.rel_SNPs]] - self.norm_const[self.dp[i+1:, self.rel_SNPs]]

            dist.append(np.nansum(logl / norm, axis=1) / valid.sum(axis=1))
        return np.concatenate(dist)


    def get_profile(self, cl):
        cells = np.isin(self.assignment, cl)
        # p = self.alt[cells].sum(axis=0) / self.dp[cells].sum(axis=0)
        p = np.average(np.nan_to_num(self.VAF[cells]), weights=self.dp[cells],
            axis=0)
        dp = np.mean(self.dp[cells], axis=0).round()
        alt = (dp * p).round()
        ref = dp - alt
        return np.hstack([alt, ref, dp])


    def check_hom_match(self, scl, dcl):
        rs =  np.reshape(self.profiles[scl], (3, -1))
        rd1 = np.reshape(self.profiles[dcl[0]], (3, -1))
        rd2 = np.reshape(self.profiles[dcl[1]], (3, -1))

        valid = (rs[2] > 0) & (rd1[2] > 0) & (rd2[2] > 0)
        rs = rs[:,valid]
        rd1 = rd1[:,valid]
        rd2 = rd2[:,valid]

        ps = rs[0] / rs[2]
        pd1 = rd1[0] / rd1[2]
        pd2 = rd2[0] / rd2[2]

        cutoff = 0.05
        hom_scl = np.argwhere(ps > 1 - cutoff**2).flatten()
        wt_scl = np.argwhere(ps <= cutoff**2).flatten()
        hom_dcl1 = set(np.argwhere(pd1 > (1 - cutoff)).flatten())
        hom_dcl2 = set(np.argwhere(pd2 > (1 - cutoff)).flatten())

        hom12_id = np.array(list(hom_dcl1 & hom_dcl2))
        hom1_id = np.array(list(hom_dcl1 - hom_dcl2))
        hom2_id = np.array(list(hom_dcl2 - hom_dcl1))

        # Check if hom in singlet match
        if hom_scl.size > 0:
            hom_s = (pd1[hom_scl] > (1 - cutoff)) & (pd2[hom_scl] > (1 - cutoff))
        else:
            hom_s = np.array([], dtype=bool)
        # Check if WT in singlet match
        if wt_scl.size > 0:
            wt_s = (pd1[wt_scl] <= cutoff) & (pd2[wt_scl] <= cutoff)
        else:
            wt_s = np.array([], dtype=bool)

        # Check if both hom match
        if hom12_id.size > 0:
            hom12 = ps[hom12_id] >= cutoff**2
        else:
            hom12 = np.array([], dtype=bool)
        # Check if hom on cl1 match
        if hom1_id.size > 0:
            hom1 = (ps[hom1_id] > cutoff**2) & (ps[hom1_id] < 1 - cutoff**2)
        else:
            hom1 = np.array([], dtype=bool)
        # Check if hom on cl2 match
        if hom2_id.size > 0:
            hom2 = (ps[hom2_id] > cutoff**2) & (ps[hom2_id] < 1 - cutoff**2)
        else:
            hom2 = np.array([], dtype=bool)

        hom_match = np.concatenate([hom_s, wt_s, hom12, hom1, hom2])

        return (~hom_match).sum()


    @staticmethod
    def get_cmap():
        cmap = LinearSegmentedColormap.from_list('my_gradient', (
            (0.000, (1.000, 1.000, 1.000)),
            (0.500, (1.000, 0.616, 0.000)),
            (1.000, (1.000, 0.000, 0.000)))
        )
        return cmap


    def get_hm_data(self):
        df =  np.nan_to_num(self.VAF, nan=-1)
        mask = np.zeros(self.VAF.shape, dtype=bool)
        mask[self.dp == 0] = True
        return df, mask


    @staticmethod
    def get_cm_specifics():
        return {'vmin': 0, 'vmax': 1,
            'cbar_kws': {'ticks': [0, 1], 'label': 'VAF'}}


    @staticmethod
    def apply_cm_specifics(cm):
        pass


    def safe_profiles(self, output):
        # Safe SNV profiles to identy patients
        cl_map, _ = self.get_cl_map()
        VAF_df = self.get_VAF_profile(self.sgt_ids)
        VAF_df.index = [f'{cl_map[i]} ({COLORS_STR[int(cl_map[i])]})' \
            for i in self.sgt_ids]
        VAF_df.round(2).to_csv(f'{output}.profiles.tsv', sep='\t')


    def get_VAF_profile(self, cl):
        reads_all = np.reshape(self.profiles[cl], (cl.size, 3, self.SNPs.size))
        VAF_raw = []
        for reads_cl in reads_all:
            VAF_raw.append(reads_cl[0] / reads_cl[2])
        return pd.DataFrame(VAF_raw, index=cl, columns=self.SNPs)


    def print_summary(self):
        cl_map, _ = self.get_cl_map()
        GTs = {'WT': (0, 0.35), 'HET': (0.35, 0.95), 'HOM': (0.95, 1)}
        for cl_id, cl_size in zip(*np.unique(self.assignment, return_counts=True)):
            cl_name = cl_map[cl_id]
            if '+' in cl_name:
                cl1, cl2 = map(int, cl_name.split('+'))
                cl_color = f'{COLORS_STR[cl1]}+{COLORS_STR[cl2]}'
            else:
                cl_color = COLORS_STR[int(cl_name)]
            print(f'Cluster {cl_name} ({cl_color}): {cl_size: >4} cells ' \
                f'({cl_size / self.cells.size * 100: >2.0f}%)')

            VAF_cl = self.VAF[self.assignment == cl_id]
            VAF_cl_called = (VAF_cl >= EPSILON).sum(axis=0)
            for geno, (geno_min, geno_max) in GTs.items():
                VAF_cl_geno = (VAF_cl >= geno_min) & (VAF_cl < geno_max)
                avg = VAF_cl_geno.sum(axis=1).mean()
                clonal = (VAF_cl_geno.sum(axis=0) > VAF_cl_called * 0.95).sum()
                print(f'\tGT: {geno} - Avg./cell {avg: >4.1f}, # 95% clonal: {clonal}')




# ------------------------------------------------------------------------------


class demoTape_gt(demoTape):
    def __init__(self, in_file, cl_no, distance='manhatten'):
        super().__init__(in_file, cl_no)
        self.prof_len = self.SNPs.size
        self.distance = distance
        if self.distance == 'manhatten':
            self.metric = self.dist_L1
        else:
            self.metric = self.dist_L2


    def __str__(self):
        return str(super()) + f'Distance metric: {self.distance}'


    def get_data(self, df):
        # Get genotype data
        self.gt = df.applymap(lambda x: float(x.split(':')[2])).values.T
        self.gt[self.gt == 3] = np.nan


    @staticmethod
    def dist_L2(c1, c2):
        valid = ~np.isnan(c1) & ~np.isnan(c2)
        return euclidean(c1[valid], c2[valid]) / np.sqrt(np.sum(2**2 * valid.sum()))


    @staticmethod
    def dist_L1(c1, c2):
        valid = ~np.isnan(c1) & ~np.isnan(c2)
        return cityblock(c1[valid], c2[valid]) / (2 * valid.sum())


    def get_pairwise_dists(self):
        dist = []
        for i in np.arange(self.cells.size - 1):
            valid = np.isfinite(self.gt[i]) & np.isfinite(self.gt[i+1:])
            manh = np.nansum(np.abs(self.gt[i] - self.gt[i+1:]), axis=1)
            norm = 2 * valid.sum(axis=1)
            dist.append(manh / norm)
        return np.concatenate(dist)


    def get_profile(self, cl):
        cells = np.isin(self.assignment, cl)
        return np.nanmean(self.gt[cells], axis=0)


    def check_hom_match(self, scl, dcl):
        hom_scl = np.argwhere(self.profiles[scl] > 1.5).flatten()
        wt_scl = np.argwhere(self.profiles[scl] <= 0.5).flatten()
        hom_dcl1 = set(np.argwhere(self.profiles[dcl[0]] > 1.5).flatten())
        hom_dcl2 = set(np.argwhere(self.profiles[dcl[1]] > 1.5).flatten())

        hom12_id = np.array(list(hom_dcl1 & hom_dcl2))
        hom1_id = np.array(list(hom_dcl1 - hom_dcl2))
        hom2_id = np.array(list(hom_dcl2 - hom_dcl1))

        # Check if hom in singlet match
        if hom_scl.size > 0:
            hom_s = (self.profiles[dcl[0]][hom_scl] > 1.5) \
                & (self.profiles[dcl[1]][hom_scl] > 1.5)
        else:
            hom_s = np.array([], dtype=bool)
        # Check if WT in singlet match
        if wt_scl.size > 0:
            wt_s =  (self.profiles[dcl[0]][wt_scl] <= 1.5) \
                & (self.profiles[dcl[1]][wt_scl] <= 1.5)
        else:
            wt_s = np.array([], dtype=bool)
        # Check if both hom match
        if hom12_id.size > 0:
            hom12 = self.profiles[scl][hom12_id] > 1.5
        else:
            hom12 = np.array([], dtype=bool)
        # Check if hom on cl1 match
        if hom1_id.size > 0:
            hom1 = (self.profiles[scl][hom1_id] >= 0.5) \
                & (self.profiles[scl][hom1_id] < 1.5)
        else:
            hom1 = np.array([], dtype=bool)
        # Check if hom on cl2 match
        if hom2_id.size > 0:
            hom2 = (self.profiles[scl][hom2_id] >= 0.5) \
                & (self.profiles[scl][hom2_id] < 1.5)
        else:
            hom2 = np.array([], dtype=bool)

        hom_match = np.concatenate([hom_s, wt_s, hom12, hom1, hom2])

        return (~hom_match).sum()


    @staticmethod
    def get_cmap():
        myColors = ('#5B566C', '#fed976', '#fc4e2a', '#800026')
        return LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))


    def get_hm_data(self):
        df = np.nan_to_num(self.gt, nan=-1)
        mask = np.zeros(self.gt.shape, dtype=bool)
        mask[np.isnan(self.gt)] = True
        return df, mask


    @staticmethod
    def get_cm_specifics():
        return {'vmin': -1, 'vmax': 2,
            'cbar_kws': {'ticks': [-0.65, 0.15, 0.85, 1.65], 'shrink': 0.5,
                'drawedges': True}}


    @staticmethod
    def apply_cm_specifics(cm):
        colorbar = cm.ax_heatmap.collections[0].colorbar
        colorbar.set_ticklabels([r' $-$', '0|0', '0|1', '1|1'])


    def safe_profiles(self, output):
        # Safe SNV profiles to identy patients
        cl_map, _ = self.get_cl_map()
        idx = [f'{cl_map[i]} ({COLORS_STR[int(cl_map[i])]})' for i in self.sgt_ids]
        pd.DataFrame(self.profiles[self.sgt_ids], index=idx, columns=self.SNPs) \
            .round(2) \
            .to_csv(f'{output}.profiles.tsv', sep='\t')



    def print_summary(self):
        cl_map, assignment_str = self.get_cl_map()
        # Print cluster summaries to stdout
        for cl_id, cl_size in zip(*np.unique(self.assignment, return_counts=True)):
            cl_name = cl_map[cl_id]
            if '+' in cl_name:
                cl1, cl2 = map(int, cl_name.split('+'))
                cl_color = f'{COLORS_STR[cl1]}+{COLORS_STR[cl2]}'
            else:
                cl_color = COLORS_STR[int(cl_name)]
            print(f'Cluster {cl_name} ({cl_color}): {cl_size: >4} cells ' \
                f'({cl_size / self.cells.size * 100: >2.0f}%)')

            gt_cl = self.gt[self.assignment == cl_id]
            gt_cl_called = np.isfinite(gt_cl).sum(axis=0)
            for geno in [0, 1, 2]:
                gt_cl_geno = gt_cl == geno
                print(f'\tGT: {geno} - ' \
                    f'Avg./cell {gt_cl_geno.sum(axis=1).mean(): >4.1f}, '
                    f'# 95% clonal: {(gt_cl_geno.sum(axis=0) > gt_cl_called * 0.95).sum()}')


# ------------------------------------------------------------------------------


def main(args):
    in_files = []
    if len(args.input) == 1 and os.path.isdir(args.input[0]):
        for file in os.listdir(args.input[0]):
            if file.endswith('_variants.csv'):
                in_files.append(os.path.join(args.input[0], file))
    else:
        in_files = args.input

    for in_file in in_files:
        if args.metric == 'reads':
            dt = demoTape_reads(in_file, args.clusters)
        else:
            dt = demoTape_gt(in_file, args.clusters, args.metric)

        dt.demultiplex()
        if args.output:
            args.output = os.path.splitext(args.output)[0]
        else:
            args.output = os.path.splitext(in_file)[0]

        dt.safe_results(args.output)
        if args.output_plot:
            dt.plot_heatmap(f'{args.output}.heatmap.{FILE_EXT}')
        if args.show_plot:
            dt.plot_heatmap()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input _variants.csv file(s).')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output base name: <DIR>/<> .')
    parser.add_argument('-n', '--clusters', type=int, default=2,
        help='Number of clusters to define. Default = 2.')
    parser.add_argument('-m', '--metric', type=str, default='reads',
        choices=['reads', 'manhatten', 'euclidean'],
        help='Which data/distance to use for demultiplexing. Default = reads.')

    plotting = parser.add_argument_group('plotting')
    plotting.add_argument('-op', '--output_plot', action='store_true',
        help='Output file for heatmap with dendrogram to "<INPUT>.hm.png".')
    plotting.add_argument('-sp', '--show_plot', action='store_true',
        help='Show heatmap with dendrogram at stdout.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)