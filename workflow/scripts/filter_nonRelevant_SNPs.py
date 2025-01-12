#!/usr/bin/env python3

import argparse
import os

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.stats import kstest, median_abs_deviation, ttest_ind
import seaborn as sns

EPSILON = np.finfo(np.float64).resolution

FILE_EXT = '.png'
FONTSIZE = 14
DPI = 300
sns.set_style('whitegrid') #darkgrid, whitegrid, dark, white, ticks
sns.set_context('paper',
    rc={'lines.linewidth': 1,
        'axes.axisbelow': True,
        'font.size': FONTSIZE,
        'axes.labelsize': 'medium',
        'axes.titlesize': 'large',
        'xtick.labelsize': 'medium',
        'xtick.major.size': 1,
        'xtick.major.width': 0.5,
        'xtick.bottom': False,
        'ytick.labelsize': 'medium',
        'legend.fontsize': 'medium',
        'legend.title_fontsize': 'large',
        'axes.labelticksize': 10
})


def plot_symmetry(vaf, vaf_ref, title=None):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(np.linspace(0, 1, vaf.size), vaf, alpha=0.75, label='VAF')
    ax.plot(np.linspace(0, 1, vaf_ref.size), vaf_ref, alpha=0.75,
        label='1 - reverse(VAF)')
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, alpha=0.25, color='grey', ls='--')
    ax.set_xlabel('rel. barcode rank')
    ax.set_ylabel('VAF')
    ax.set_title(title)
    ax.legend()
    plt.show()


def load_variants(in_file):
    return pd.read_csv(in_file, dtype={'CHR': str})


def load_cell_annotation(in_file):
    df = pd.read_csv(in_file, index_col=0)
    if 'assignmnet' in df.columns:
        df.rename({'assignmnet': 'assignment'}, axis=1, inplace=True)
    return df


def test_for_relevant_SNPs(df, df_ca, dbt_cells):
    snp_names = df.apply(lambda x: f'{x["CHR"]}:{x["POS"]} {x["REF"]}>{x["ALT"]}',
        axis=1)
    vals = df.iloc[:, 7:].drop(dbt_cells, axis=1)
    ref = vals.map(lambda x: int(x.split(':')[0]))
    alt = vals.map(lambda x: int(x.split(':')[1]))
    gt = vals.map(lambda x: int(x.split(':')[2]))
    dp = alt + ref
    VAF = (alt + EPSILON) / (dp + 2 * EPSILON)

    # nans = vals.map(lambda x: x.split(':')[2] == '3')
    nans = dp < 2
    VAF = VAF.mask(nans)

    cols = ['>0.95 VAF', '<0.05 VAF', 'VAF std', 'symmetry kstest q-value']
    if df_ca.size > 0:
        h_idx = df_ca[df_ca['assignment'] == 'healthy'].index
        t_idx = df_ca[df_ca['assignment'] == 'tumor'].index
        
        df_h = VAF[h_idx]
        df_t = VAF[t_idx]

        dp_h = dp[h_idx]
        dp_t = dp[t_idx]

        cols.extend(['groups kstest q-value', 'groups kstest statistic',
            'mean VAF tumor', 'mean VAF healthy', 'mean VAF FC',
            'std VAF tumor', 'std VAF healthy', 'std VAF FC', 
            'median VAF tumor', 'median VAF healthy', 'median VAF FC',
            'MAD VAF tumor', 'MAD VAF healthy', 'MAD VAF FC', 
            'perc do tumor', 'groups do ttest q-value'])

    no_nans = (~nans).sum(axis=1)
    VAF_hom = (VAF > 0.95).sum(axis=1) / no_nans
    VAF_wt = (VAF < 0.05).sum(axis=1) / no_nans
    VAF_std = VAF.std(axis=1)
    
    data = []
    # Iterate over SNPs
    for snp, snp_data in VAF.iterrows():
        vals = snp_data[dp.loc[snp] >= 10].values
        # Add correlation between sorted VAF distribution and 
        #    1 - reverse sorted VAF distribution: symmetric/normal distributed 
        #        --> likely germline + ADO

        # Possible allelic dropout: only 1 allele amplified or too low depth
        pos_ado = (vals <= 0.001) | (vals >= 0.999)

        low_vafs = vals[(vals < 0.5) & ~pos_ado]
        high_vafs = 1 - vals[(vals > 0.5) & ~pos_ado]

        # VAFs are very unevenly distributed anyhow: q-Value = 0
        if low_vafs.size < 6 or high_vafs.size < 6:
            ks_sym_qval = np.nan
        # CHANGED: additional filter added after demoTape submission
        # Only 50% of data with intermediate VAFs: low q-Value
        elif pos_ado.mean() > 0.5:
            ks_sym_qval = np.nan
        else:
            ks_sym_pval = kstest(low_vafs, high_vafs, nan_policy='omit').pvalue
            ks_sym_qval = min(1, ks_sym_pval * snp_names.size)

        # plot_symmetry(low_vafs, high_vafs, snp_names[snp])

        snp_res = [VAF_hom[snp], VAF_wt[snp], VAF_std[snp], ks_sym_qval]

        if df_ca.size > 0:
            ks_grp = kstest(df_h.loc[snp], df_t.loc[snp], nan_policy='omit')
            ks_grp_qval = min(1, ks_grp.pvalue * VAF.shape[1])

            t_mean = df_t.loc[snp].mean()
            h_mean = df_h.loc[snp].mean()
            t_std = df_t.loc[snp].std()
            h_std = df_h.loc[snp].std()
            t_median = df_t.loc[snp].median()
            h_median = df_h.loc[snp].median()
            t_mad = median_abs_deviation(df_t.loc[snp], nan_policy='omit')
            h_mad = median_abs_deviation(df_h.loc[snp], nan_policy='omit')

            t_perc_nan = (dp_t.loc[snp] < 2).mean()
            ttest_nan = ttest_ind(
                (dp_t.loc[snp] < 2).astype(int), (dp_h.loc[snp] < 2).astype(int))
            ttest_nan_qval = min(1, ttest_nan.pvalue * VAF.shape[1])

            snp_res.extend([
                ks_grp_qval, ks_grp.statistic,
                t_mean, h_mean, (t_mean + EPSILON) / (h_mean + EPSILON),
                t_std, h_std, (t_std + EPSILON) / (h_std + EPSILON),
                t_median, h_median, (t_median + EPSILON) / (h_median + EPSILON),
                t_mad, h_mad, (t_mad + EPSILON) / (h_mad + EPSILON), 
                t_perc_nan, ttest_nan_qval])

        data.append(snp_res)
    df_out = pd.DataFrame(data, index=VAF.index, columns=cols)

    return df_out


def check_against_dbsnp(in_file, df_rel, clinVar=False):
    df = pd.read_csv(in_file, sep='\t')
    df['CHROM'] = df['CHROM'].str.replace('chr', '')
    df.set_index(['CHROM', 'POS'], inplace=True)
    df = df.sort_index()

    if clinVar:
        col_name = 'clinVar het freq'
    else:
        col_name = 'dbsnp het freq'

    df_rel[col_name] = 0.0
    for row, data in df_rel.iterrows():
        try:
            dbsnp_entry = df.loc[row]
        except (KeyError, TypeError):
            pass
        else:
            if dbsnp_entry.loc[row, 'ALT'] == data['ALT']:
                df_rel.at[row, col_name] = dbsnp_entry.loc[row, 'EXPECTEDHET']


def check_fc(x, cutoff=0.1):
    return (1 - cutoff < x) & (x < 1 + cutoff)


def main(args):
    if not args.output:
        args.output = args.variants
    base_dir = os.path.dirname(args.output)
    prefix = os.path.basename(args.output).split('.filtered_variants')[0]

    df_v = load_variants(args.variants)
    df_g = load_variants(args.genotypes)

    if args.cell_annotation:
        df_ca = load_cell_annotation(args.cell_annotation)
        if 'doublets' in df_ca['assignment'].unique():
            dbt_cells = df_ca[df_ca['assignment'] == 'doublets'].index
        else:
            dbt_cells = []
    else:
        df_ca = pd.DataFrame([])
        dbt_cells = []

    df_rel = test_for_relevant_SNPs(df_v, df_ca, dbt_cells)
    # Add CHR, POS, REF, ALT columns
    id_cols = ['CHR', 'POS', 'REF', 'ALT']
    df_rel[id_cols] = df_v.loc[:,id_cols]
    df_rel.set_index(['CHR', 'POS'], inplace=True)
    # Reorder columns such that REF and ALT are first and second
    df_rel = df_rel[['REF', 'ALT'] + df_rel.columns[:-2].tolist()]

    if args.dbsnp:
        check_against_dbsnp(args.dbsnp, df_rel)
    if args.clinVar:
        check_against_dbsnp(args.clinVar, df_rel, clinVar=True)
    df_rel['filter'] = ''

    # Select only relevant SNPs
    if args.cell_annotation:
        # Make sure there is no group difference in nans, as CN 0 could cause this
        hom_sides = (df_rel['>0.95 VAF'] >= 0.99) \
            & ((df_rel['perc do tumor'] <= 0.05) \
                | (df_rel['groups do ttest q-value'] > 0.05))
        df_rel.loc[hom_sides, 'filter'] += ';all hom'
    else:
        df_rel.loc[(df_rel['>0.95 VAF'] >= 0.99), 'filter'] += ';all hom'
    df_rel.loc[(df_rel['<0.05 VAF'] >= 0.99), 'filter'] += ';all wt'
    df_rel.loc[df_rel['symmetry kstest q-value'] > args.alpha, 'filter'] \
        += ';VAFsymmetry'

    if args.dbsnp:
        df_rel.loc[df_rel['dbsnp het freq'] > 0.01, 'filter'] += ';dbsnp'
    if args.clinVar:
        df_rel.loc[df_rel['clinVar het freq'] > 0, 'filter'] += ';clinVar'
    if args.cell_annotation:
        no_diff = df_rel['groups kstest q-value'] > args.alpha
        no_effect = check_fc(df_rel['mean VAF FC']) & check_fc(df_rel['std VAF FC']) \
            & check_fc(df_rel['median VAF FC']) & check_fc(df_rel['MAD VAF FC']) 
        df_rel.loc[no_diff | no_effect, 'filter'] \
            += ';no grp diff'
    
        germline_snp = df_rel['filter'].str.contains('dbsnp') \
            & (df_rel['mean VAF healthy'] > 0.33) \
            & check_fc(df_rel['mean VAF FC']) \
            & (df_rel['median VAF healthy'] > 0.33) \
            & check_fc(df_rel['MAD VAF FC'])
        df_rel.loc[germline_snp, 'filter'] += ';germline'

    df_rel['filter'] = df_rel['filter'].str.lstrip(';')

    keep_idx = np.argwhere(
            ((df_rel['filter'] == '') \
                | (df_rel['filter'].str.contains('clinVar')) \
                | ((df_rel['filter'] == 'dbsnp') & (df_rel['VAF std'] >= 0.05))
            ).values) \
        .ravel()

    var_out = os.path.join(base_dir, f'{prefix}.relevant.filtered_variants.csv')
    df_v.iloc[keep_idx].to_csv(var_out, index=False)
    gt_out = os.path.join(base_dir, f'{prefix}.relevant.filtered_variants_gt.csv')
    df_g.iloc[keep_idx].to_csv(gt_out, index=False)

    stats_file = os.path.join(base_dir, f'{prefix}.relevant.statistics.csv')
    print(f'Writing statistics to: {stats_file}')
    df_rel.to_csv(stats_file)

    # Create data for SNP vulcano style plot:
    #   y axis = p-values: Kolmogorov-Smirnov-test
    #   x axis = effect size: mean|median diff | ks test statistic
    if df_ca.size > 0:
        out_vulcano =  os.path.join(base_dir, f'{prefix}.SNP_vulcano.png')
        plot_vulcano(df_rel, out_vulcano)


def plot_vulcano(df, out_file):
    nrows = 1
    ncols = 3
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
        figsize=(5 * ncols, 5 * nrows))
    axes = np.reshape(axes, (nrows, ncols))

    y = -1 * np.log10(np.clip(df['groups kstest q-value'], EPSILON, None))
    y_sig = -1 * np.log(0.05)
    for col, x_col in enumerate(
                ['mean VAF FC', 'median VAF FC', 'groups kstest statistic']):
        ax = axes[0][col]
        x = df[x_col]
        if x_col == 'groups kstest statistic':
            down_diff = pd.Series([False] * x.shape[0], index=x.index)
            up_diff = (y > y_sig) & (x > 0.05)
        else:
            down_diff = (y > y_sig) & (x < 0.9)
            up_diff = (y > y_sig) & (x > 1.1)
            x = -1 * np.log10(x)

        rel = pd.Series(['normal'] * df.shape[0], index=df.index, name='relevance')
        rel.loc[down_diff] = 'down'
        rel.loc[up_diff] = 'up'
        col_map = {'normal': '#808080', 'up': '#AA3939', 'down': '#2E4172'}
        
        sp = sns.scatterplot(
            ax=ax,
            y=y,
            x=x,
            hue=rel,
            palette=col_map
        )
        ax.axhline(y=y_sig, ls='--', c='black')

        sp.set_ylabel('-log10(q-value)')
        if x_col == 'groups kstest statistic':
            sp.set_xlabel(x_col)
            sp.set_xlim([0, 1])
        else:
            ax.axvline(x=-1*np.log10(0.9), ls='--', c='black')
            ax.axvline(x=-1*np.log10(1.1), ls='--', c='black')
            sp.set_xlabel(f'-log10({x_col})')
            x_max = x.abs().max() * 1.05
            sp.set_xlim([-x_max, x_max])

        sns.move_legend(sp, 'lower right') 

    if out_file:
        fig.savefig(out_file, dpi=DPI)
    else:
        plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--variants', type=str,
        help='Variant file from mosaic preprocessing.'),
    parser.add_argument('-g', '--genotypes', type=str,
        help='Genotype file from mosaic preprocessing.'),
    parser.add_argument('-ca', '--cell_annotation', type=str, default='',
        help='Cell order (barcodes) in csv format.')
    parser.add_argument('-d', '--dbsnp', type=str, default='',
        help='DBSNP .txt file (from tapestri pipeline).')
    parser.add_argument('-c', '--clinVar', type=str, default='',
        help='clinVar .txt file (from tapestri pipeline).')
    parser.add_argument('-o', '--output', type=str, default='',
        help='"Base" file for outputs. Default: <VARIANTS_DIR>/<PREFIX>')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
        help='Significance level. Default: 0.05')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)