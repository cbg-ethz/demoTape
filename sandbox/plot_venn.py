#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2, venn3_circles, venn2_circles
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_context('talk')

import venn


SID_MAP = {'SID3841': 'S1', 'SID3867': 'S2','SID5962': 'S3',
    'mixed_SID3867_SID5962_SID3841': 'MS1',
    'SID3841 (multiplexed)': 'S1 (multiplexed)', 'SID3841 (cluster)': 'S1 (cluster)',
    'SID3867 (multiplexed)': 'S2 (multiplexed)', 'SID3867 (cluster)': 'S2 (cluster)',
    'SID5962 (multiplexed)': 'S3 (multiplexed)', 'SID5962 (cluster)': 'S3 (cluster)',
    }
CL_MAP = {
    '0': 'SID5962',
    '1': 'SID3841',
    '2': 'SID3867'
}
FIGSIZE = (10, 10)
FONTSIZE = 40
DPI = 300


CHR_MAP = dict({f'{i}': i for i in range(1,23,1)}, **{'X': 23, 'Y': 24})

default_colors = [
    # r, g, b, a
    [92, 192, 98, 0.5],
    [90, 155, 212, 0.5],
    [246, 236, 86, 0.6],
    [241, 90, 96, 0.4],
    [255, 117, 0, 0.3],
    [82, 82, 190, 0.2],
]
default_colors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in default_colors
]


def get_sets_advanced(df, cutoff):
    Abc = set([])
    aBc = set([])
    ABc = set([])
    abC = set([])
    AbC = set([])
    aBC = set([])
    ABC = set([])

    for snv_id, data in df.groupby(level=0):
        mt_freq = data.sum()
        mt_clonal = (mt_freq > cutoff).astype(int)
        ab_clonal = ((mt_freq <= max(0, (1 - cutoff))) | (mt_freq.isnull())).astype(int)

        clonal = (data > cutoff).sum()
        subclonal = ((data <= cutoff) & (data > (1 - cutoff))).sum()
        absent = ((data <= max(0, (1 - cutoff))) | (data.isnull())).sum()
        # No freq about cutoff in any sample
        if not clonal.any():
            continue

        cl_names = np.where(clonal)[0]
        subcl_names = np.where(subclonal)[0]
        mt_names = np.concatenate([cl_names, subcl_names])
        ab_names = np.where(absent == 2)[0]

        if cl_names.size == 1:
            # Only clonal in 1 sample
            if ab_names.size == 2:
                if cl_names[0] == 0:
                    Abc.add(snv_id)
                elif cl_names[0] == 1:
                    aBc.add(snv_id)
                else:
                    abC.add(snv_id)
            # Also subclonal in other samples
            else:
                if mt_names.size == 3:
                    ABC.add(snv_id)
                elif mt_names.size == 2:
                    if 0 in mt_names and 1 in mt_names:
                        ABc.add(snv_id)
                    elif 0 in mt_names and 2 in mt_names:
                        ABc.add(snv_id)
                    else:
                        aBC.add(snv_id)
                else:
                    print(f'Unclear:\n{data}\n')

        elif cl_names.size == 2:
            if ab_names.size == 1:
                if 0 in cl_names and 1 in cl_names:
                    ABc.add(snv_id)
                elif 0 in cl_names and 2 in cl_names:
                    ABc.add(snv_id)
                else:
                    aBC.add(snv_id)
            else:
                if mt_names.size == 3:
                    ABC.add(snv_id)
                else:
                    print(f'Unclear:\n{data}\n')
        else:
            ABC.add(snv_id)
    return ([Abc, aBc, ABc, abC, AbC, aBC, ABC])


def bit_to_str(bit_str, labels):
    bit_map = {}
    for i in bit_str:
        new_str = ''
        for j, k in enumerate(i):
            if k == '1':
                new_str += labels[j] + '+'
        bit_map[i] = new_str[:-1]
    return bit_map


def plot_venn(freq_df, out_file, cutoff, specific):
    names = freq_df.columns.values
    names = [SID_MAP[i] if i in SID_MAP else i for i in names]

    SNPs = []
    for i in freq_df.columns:
        clonal = freq_df[freq_df[i] > cutoff] \
            .index.get_level_values(0).unique().values
        if args.specific:
            rel_cols = freq_df.loc[:, freq_df.columns != i]
            spec = freq_df[rel_cols.max(axis=1) < (1 - cutoff)] \
                .index.get_level_values(0).unique().values
        else:
            spec = clonal
        SNPs.append(list(set(clonal) & set(spec)))

    labels, cnt_map = venn.get_labels(SNPs)
    sampel_map = bit_to_str(labels.keys(), names)
    cnt_map = {sampel_map[i]: j for i,j in cnt_map.items()}

    if len(names) == 3:
        counts = np.array([int(labels['100']), int(labels['010']),
            int(labels['110']), int(labels['001']), int(labels['101']),
            int(labels['011']), int(labels['111'])])
        plt_fct = venn3
        circle_fct = venn3_circles
    elif len(names) == 2:
        counts = np.array([int(labels['10']), int(labels['01']), int(labels['11'])])
        plt_fct = venn2
        circle_fct = venn2_circles
    elif len(names) == 4:
        fig, ax = venn.venn4(labels, names=names, fontsize=FONTSIZE, dpi=DPI,
            figsize=FIGSIZE)
        ax.get_legend().remove()
        print('Saving Venn-plots to: {}'.format(out_file))
        fig.tight_layout()
        fig.savefig(out_file, dpi=300)
        return cnt_map
    else:
        raise RuntimeError(f'Unknown number of circles to draw: {len(names)}')

    counts_norm = np.clip(counts / counts.sum(), 0.025, 1).round(2)
    if len(names) == 3:
        counts_norm = counts_norm + np.arange(0.0001, 0.00071, 0.0001)
    else:
        counts_norm = counts_norm + np.arange(0.0001, 0.00031, 0.0001)

    no_formatter = {}
    for i, j in enumerate(counts_norm):
        no_formatter[j] = counts[i]

    def formatter(x):
        return no_formatter[x]
    # ------------------------------

    fig = plt.figure(figsize=FIGSIZE)
    v = plt_fct(
        subsets=(counts_norm),
        set_labels=names,
        set_colors=default_colors[:len(names)],
        alpha=0.5,
        normalize_to=counts_norm.sum(),
        subset_label_formatter=formatter
    )
    c = circle_fct(
        subsets = (counts_norm),
        normalize_to=counts_norm.sum()
    )

    for text in v.set_labels:
        if text:
            text.set_fontsize(FONTSIZE + 10)

    for text in v.subset_labels:
        if text:
            text.set_fontsize(FONTSIZE)


    try:
        fig.tight_layout()
    except AttributeError:
        pass

    print('Saving Venn-plots to: {}'.format(out_file))
    fig.savefig(out_file, dpi=300)
    plt.close()
    return cnt_map


def main(args):
    if args.input_variant:
        freq_df = get_SNVs_variant(args.input_variant, args.cluster_assignment)
    elif args.input_profiles:
        freq_df =  get_SNVs_profiles(args.input_profiles, ids, names, args.cutoff)
    else:
        raise IOError('No input file(s) given (profile file | loom files)')

    out_str_raw = "+".join(freq_df.columns.values)
    if args.specific:
        freq_df.dropna(inplace=True)
        out_str_raw += '_spec'

    out_str_cl = f'{out_str_raw}_{float(args.cutoff)*100:0>2.0f}clonal'


    if not args.out_file:
        if args.input_variant:
            out_dir = os.path.dirname(args.input_variant[0])
        elif args.input_profiles:
            out_dir = os.path.dirname(args.input_profiles[0])

        args.out_file = os.path.join(out_dir, f'{out_str_cl}.png')
    else:
        out_dir = os.path.dirname(args.out_file)

    cnt_map = plot_venn(freq_df, args.out_file, args.cutoff, args.specific)
    data_file = os.path.join(out_dir, f'{out_str_cl}.tsv')
    safe_table(cnt_map, data_file)
    plot_freq_profiles(freq_df, os.path.join(out_dir, f'{out_str_raw}_freq.png'))


def plot_freq_profiles(df, out_file):
    cmap = plt.get_cmap('YlOrRd', 100)
    cmap.set_under('#EAEAEA')

    data_plot = []
    labels = []
    SNPs = np.unique(df.index.get_level_values(0).values)

    for zyg in ['het', 'hom']:
        for col in sorted(df.columns):
            if col in SID_MAP:
                labels.append(f'{SID_MAP[col]} {zyg}. ')
            else:
                labels.append(f'{col} {zyg}.')
            data_plot.append(df.xs(zyg, level=1).loc[SNPs, col].fillna(-1).values)
        labels.append('')
        data_plot.append(np.full(SNPs.size, np.nan))

    df_plot = pd.DataFrame(data_plot[:-1], index=labels[:-1], columns=SNPs)
    SNPs_sorted = sorted(SNPs,
        key=lambda x: (CHR_MAP[x.split('_')[0]], int(x.split('_')[1])))
    df_plot = df_plot[SNPs_sorted]
    df_plot.columns = ['{}:{} {}>{}'.format(*i.split('_')) for i in SNPs_sorted]

    fig, ax = plt.subplots(figsize=(SNPs.size // 2, 10))

    hm = sns.heatmap(
        df_plot,
        annot=False,
        square=True,
        cmap=cmap,
        vmin=0,
        vmax=1,
        cbar_kws={'label': 'Frequency', 'shrink': 0.5, 'ticks': [0, 0.5, 1]},
        linewidths=0,
        ax=ax
    )
    ax.set_xlabel('SNPs')
    ax.set_ylabel('')
    ax.yaxis.get_ticklabels()[3].set_visible(False)

    fig.tight_layout()
    if out_file:
        print(f'Frequency heatmap written to: {out_file}')
        fig.savefig(out_file, dpi=300)
    else:
        plt.show()



def safe_table(cnt_map, out_file):
    print('Saving Venn-data to: {}'.format(out_file))
    with open(out_file, 'w') as f:
        f.write('Sample\tSNVs')
        for sample, SNVs in cnt_map.items():
            SNV_str = ';'.join(sorted(list(SNVs)))
            f.write(f'\n{sample}\t{SNV_str}')


def get_SNVs_variant(in_files, clusters=''):
    if clusters:
        # If multiplexed and single SNP profiles, multiplexed has to be first!
        freq_df = get_SNVs_multiplex(in_files[0], clusters)
        if len(in_files) > 1:
            freq_df2 = get_SNVs_single(in_files[1:])
            # If second file is single sample, take only this one
            if freq_df2.columns[0] in SID_MAP:
                freq_df2.rename({freq_df2.columns[1]: freq_df2.columns[0] + ' (cluster)'}, axis=1, inplace=True)
                freq_df = freq_df[freq_df2.columns[0]]
                mltpl_col = freq_df.name + ' (multiplexed)'
                freq_df.name = mltpl_col
            freq_df = pd.DataFrame(freq_df).merge(freq_df2,
                left_index=True, right_index=True, how='outer')
            freq_df = freq_df[list(freq_df2.columns.values) + [mltpl_col]]
    else:
        freq_df = get_SNVs_single(in_files)
    return freq_df


def get_SNVs_single(in_files):
    for in_file in in_files:
        df = pd.read_csv(in_file, dtype=str)
        df.dropna(inplace=True)
        df['idx'] = df['CHR'] + '_' + df['POS'] + '_' + df['REF'] + '_' + df['ALT']
        df.set_index('idx', inplace=True)

        gt = df.iloc[:,7:].applymap(lambda x: int(x.split(':')[-1])).values

        het_freq = (gt == 1).sum(axis=1) / (gt != 3).sum(axis=1)
        hom_freq = (gt == 2).sum(axis=1) / (gt != 3).sum(axis=1)

        path_split = in_file.split(os.sep)
        if in_file.endswith('_variants.csv'):
            name = path_split[-1].replace('.cells_variants.csv', '')
        else:
            if re.match('\d\d-\d\d', path_split[-2]):
                file_name = path_split[-1]
                if in_file.endswith('.filtered_variants.csv'):
                    mix = '{}:{}'.format(*path_split[-2].split('-'))
                    name = file_name.replace('.filtered_variants.csv', f' {mix}')
                else:
                    name = file_name.split('_')[0]
            else:
                name = path_split[-2]

        idx_tuples = list(zip(
            np.repeat([df.index.values], 2, axis=0).flatten(),
            np.repeat(['het', 'hom'], df.shape[0])))
        index = pd.MultiIndex.from_tuples(idx_tuples, names=['id', 'zygosity'])
        new_freq_df = pd.DataFrame(np.concatenate([het_freq, hom_freq]),
            index=index, columns=[name])

        try:
            freq_df = freq_df.merge(new_freq_df,
                left_index=True, right_index=True, how='outer')
        except NameError:
            freq_df = new_freq_df

    return freq_df



def get_SNVs_multiplex(in_file, clusters):
    df = pd.read_csv(in_file, dtype=str)
    df['idx'] = df['CHR'] + '_' + df['POS'] + '_' + df['REF'] + '_' + df['ALT']
    df.set_index('idx', inplace=True)
    cells = df.columns[7:].values

    assignment = pd.read_csv(clusters, index_col=0, dtype=str, sep='\t')
    assign_vec = assignment.values.flatten()

    for cluster in np.unique(assign_vec):
        if '+' in cluster:
            continue

        try:
            name = CL_MAP[cluster]
        except KeyError:
            name = cluster

        cl_cells = assignment.columns[np.argwhere(assign_vec == cluster).flatten()].values
        cl_df = df[cl_cells]
        cl_gt = np.array(cl_df.applymap(lambda x: int(x.split(':')[-1])).values)

        het_freq = (cl_gt == 1).sum(axis=1) / (cl_gt != 3).sum(axis=1)
        hom_freq = (cl_gt == 2).sum(axis=1) / (cl_gt != 3).sum(axis=1)

        cl_freq = ((cl_gt == 1) | (cl_gt == 2)).sum(axis=1) / (cl_gt != 3).sum(axis=1)

        idx_tuples = list(zip(
            np.repeat([cl_df.index.values], 2, axis=0).flatten(),
            np.repeat(['het', 'hom'], cl_df.shape[0])))
        index = pd.MultiIndex.from_tuples(idx_tuples, names=['id', 'zygosity'])
        new_freq_df = pd.DataFrame(np.concatenate([het_freq,hom_freq]),
            index=index, columns=[name])

        try:
            freq_df = freq_df.merge(new_freq_df,
                left_index=True, right_index=True, how='outer')
        except NameError:
            freq_df = new_freq_df

    return freq_df


def get_SNVs_profiles(profiles_files, ids, names, cutoff=0.05):
    for profiles_file in profiles_files:
        df = pd.read_csv(profiles_file, sep='\t', index_col=0)
        # 1 File containing multiple profiles (per row)
        if df.shape[0] > 1:
            for name, data in df.iterrows():
                names.append(name)
                ids.append(set(df.columns[data > cutoff].values))

        # 1 profiles per file
        else:
            ids.append(set(df.columns[df.iloc[0] > cutoff].values))
            file_name = os.path.basename(profiles_file)
            file_name.replace('_variants.csv.profiles.tsv', '')
            names.append(file_name.split('_')[0])


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-iv', '--input_variant', type=str, nargs='+',
        help='3 different variant files to compare.')
    parser.add_argument('-ip', '--input_profiles', type=str, nargs='+',
        help='File with SNP profiles of 3 sample or 3 individual profile files.')
    parser.add_argument('-cl', '--cluster_assignment', type=str, default='',
        help='Cluster assignment for multiplex variant file.')
    parser.add_argument('-c', '--cutoff', type=float, default=0.00,
        help='Clonality cutoff: SNPs smaller than <cutoff> are removed.')
    parser.add_argument('-s', '--specific', action='store_true',
        help='Find patient specific SNPs present >= <cutoff> fraction of cells.')
    parser.add_argument('-o', '--out_file', type=str, help='Output directory')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)