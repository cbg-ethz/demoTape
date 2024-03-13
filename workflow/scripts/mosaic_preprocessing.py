#!/usr/bin/env python3

import argparse
import gc
from itertools import combinations
import os
import re
import shutil
import tempfile


import loompy
import numpy as np
import pandas as pd


CHR_ORDER = {str(i): i for i in range(1, 23)}
CHR_ORDER.update({'X':23, 'Y':24})

WHITELIST = []

VCF_HEADER ='##fileformat=VCFv4.3\n' \
    '##source=demoTape simulations\n' \
    '{contigs}' \
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Called ML genotype">\n' \
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n' \
    '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Alternative depth">\n' \
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{cells}\n'

VCF_BODY_ROW = '{chr}\t{pos}\t{chr}:{pos}\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP:AD\t{data}\n'


def plot_heatmap(df, SNPs=None, cells=None, out_file=None,
            cluster=True):
    COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#E4761A' , '#2FBF85']
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from seaborn import clustermap

    myColors = ('#EAEAEA', '#fed976', '#fc4e2a', '#800026')
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    if isinstance(cells, (list, pd.Series, np.ndarray)):
        r_colors = [[], []]
        for i in cells:
            if '+' in i:
                s1, s2 = sorted([int(j[-1]) for j in i.split('+')])
                r_colors[0].append(COLORS[s1])
                r_colors[1].append(COLORS[s2])
            else:
                r_colors[0].append(COLORS[int(i[-1])])
                r_colors[1].append('#FFFFFF')
        y_labels = cells
    else:
        r_colors = None
        y_labels = 'auto'

    if isinstance(SNPs, (list, pd.Series, np.ndarray)):
        x_labels = SNPs
    else:
        x_labels = 'auto'

    df[df == 3] = -1
    if df.shape[0] < df.shape[1]: # Wrong dimensions: transpose
        df = df.T
    cm = clustermap(
        df,
        row_cluster=cluster,
        col_cluster=cluster,
        vmin=-1, vmax=2,
        row_colors=r_colors,
        cmap=cmap,
        figsize=(25, 10),
        xticklabels=x_labels,
        yticklabels=y_labels,
        cbar_kws={'shrink': 0.5, 'drawedges': True}
    )
    cm.ax_heatmap.set_ylabel('Cells')
    if df.shape[0] > 50:
        cm.ax_heatmap.set_yticks([])
    cm.ax_heatmap.set_xlabel('SNPs')

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


def get_assignment_dict(assignment_file):
    assign = pd.read_csv(assignment_file, index_col=0, dtype='str',sep='\t')

    clusters = {}
    for cl in np.unique(assign):
        if '+' in cl:
            continue
        cl_name = f'.{cl}'
        clusters[cl_name] = assign.columns[(assign == cl).values.flatten()].values
    return clusters


def get_WL_ids(df):
    ids = df['CHR'].astype(str) + '_' + df['POS'].astype(str)
    return np.argwhere(np.isin(ids.values, WHITELIST)).flatten()


def main(args):
    if args.assignment:
        assign = get_assignment_dict(args.assignment)
    else:
        assign = {'': np.array([])}

    for cl, cells in assign.items():
        if len(args.input) == 1:
            df1, gt1, VAF = process_loom(args, cells)
        else:
            df1, gt1, VAF = multiplex_looms(args)

        # Sort dataframe by chr:pos
        df1['CHR_ORDER'] = df1['CHR'].map(CHR_ORDER)
        df1.sort_values(['CHR_ORDER', 'POS'], inplace=True)
        df1.drop('CHR_ORDER', axis=1, inplace=True)
        gt1 = gt1[df1.index.values]
        VAF = VAF[df1.index.values]
        df1.reset_index(drop=True, inplace=True)

        df2, gt2 = filter_variants(df1, gt1, VAF, args)
        df3, gt3 = filter_variants_consecutive(df2, gt2, args.proximity)
        if args.filter_same:
            df4, gt4 = postprocess_mosaic(df2, gt2)
        else:
            df4, gt4 = df3, gt3

        if len(WHITELIST) != 0:
            df4, gt4 = filter_variants_WL(df4, gt4, args)

        # Generate regions file if panel with gene name data is provided
        if args.panel:
            try:
                gene_file = [i for i in os.listdir(args.panel) \
                    if i.endswith('-submitted.bed')][0]
            except KeyError:
                raise IOError(f'No <PANEL>-submitted.bed file in {args.panel}')
                
            pa_raw = pd.read_csv(os.path.join(args.panel, gene_file),
                sep='\t', header=None, skiprows=[0])
            pa_raw[0] = pa_raw[0].str.replace('chr', '')
            regions = []
            for gene, reg_df in pa_raw.groupby(3):
                if not gene.startswith('chr'):
                    gene = gene.split(':')[0]
                regions.append(
                    [reg_df.iloc[0, 0], reg_df[1].min(), reg_df[2].max(), gene])
            pa = pd.DataFrame(regions, columns=['chr', 'start', 'stop', 'gene'])

            for i, data in df4.iterrows():
                overlap = pa[(pa['chr'] == data['CHR']) \
                    & (data['POS'] >= pa['start']) & (data['POS'] <= pa['stop'])]
                if overlap.size == 0:
                    continue
                elif overlap.shape[0] == 1:
                    df4.loc[i, 'REGION'] = overlap.iloc[0, 3]
                elif overlap['gene'].unique().size == 1:
                    df4.loc[i, 'REGION'] = overlap['gene'].unique()[0]
                elif len([i for i in overlap['gene'].unique() \
                        if not i.startswith('chr')]) == 1:
                    df4.loc[i, 'REGION'] = [i for i in overlap['gene'].unique() \
                        if not i.startswith('chr')][0]
                else:
                    import pdb; pdb.set_trace()
            map_r = {j["gene"]: f'{j["chr"]}_{j["gene"]}' for _, j in pa.iterrows()}
            for i in df4['REGION']:
                if i not in map_r:
                    chrom = df4[df4['NAME'] == i].iloc[0, 0]
                    map_r[i] = f'{chrom}_{i}'

            DP = df4.iloc[:,7:].applymap(
                lambda x: sum([i for i in map(int, x.split(':')[:2])]))
            DP['REGION'] = df4['REGION']

            DP_r = DP.groupby('REGION').sum()
            DP_r.rename(index=map_r, inplace=True)

        if args.output:
            if args.output.endswith('_variants.csv'):
                variant_file = f'{args.output[:-13]}{cl}_variants.csv'
                out_file_raw = os.path.splitext(variant_file)
            else:
                out_file_raw = os.path.splitext(args.output)
                variant_file = f'{out_file_raw[0]}{cl}{out_file_raw[1]}'
        else:
            variant_file = f'{os.path.splitext(args.input[0])[0]}{cl}_variants.csv'
            out_file_raw = os.path.splitext(variant_file)

        out_dir = os.path.dirname(variant_file)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            print(f'Creating output directory: {out_dir}')

        # Filtering germlines
        gt = gt4.astype(float)
        gt[gt == 3] = np.nan
        germline = np.nanstd(gt == 1, axis=1) < args.minStd

        loci = df4['CHR'] + ':' + df4['POS'].astype(str) + ':' + df4['REF'] \
            + '/' + df4['ALT']

        gl_file = f'{out_file_raw[0]}{cl}_germlineOrArtifact.csv'
        print(f'Writing germlines/technical artifacts to: {gl_file}')
        with open(gl_file, 'w') as f:
            f.write('\n'.join(loci[germline].values))

        gt4 = gt4[~germline]
        df4 = df4.iloc[~germline]

        print(f'Writing variant file to: {variant_file}')
        df4.to_csv(variant_file, index=False, header=True)
        if args.panel:
            region_file = variant_file.replace('_variants.csv', '_regions.csv')
            print(f'Writing region file to: {region_file}')
            DP_r.to_csv(region_file, index=True, header=False)

        if not args.full_output:
            continue
        out_file_raw = os.path.splitext(variant_file)
        print(out_file_raw)
        variant_full_file = f'{out_file_raw[0]}_full{out_file_raw[1]}'
        df1.to_csv(variant_full_file, index=False, header=True)

        save_vcf(df4.copy(deep=True), out_file_raw[0])

        cells = df4.columns[7:].values

        barcode_file = f'{out_file_raw[0]}_barcodes.csv'
        with open(barcode_file, 'w') as f:
            f.write('\n'.join(cells))

        RD = df4.iloc[:,7:].applymap(lambda x: int(x.split(':')[0])).set_index(loci)
        AD = df4.iloc[:,7:].applymap(lambda x: int(x.split(':')[1])).set_index(loci)
        DP = RD + AD

        RD_sm_file = f'{out_file_raw[0]}_RD.mtx'
        AD_sm_file = f'{out_file_raw[0]}_AD.mtx'
        with open(RD_sm_file, 'w') as f_r, open(AD_sm_file, 'w') as f_a:
            for f in [f_r, f_a]:
                f.write('%%MatrixMarket matrix coordinate real general\n' \
                    '% written by distance-demulti\n'\
                    f'{loci.size} {cells.size} {gt4.size}\n')
            for snv in range(loci.size):
                for cell in range(cells.size):
                    f_r.write(f'{snv+1} {cell+1} {RD.iloc[snv, cell]}\n')
                    f_a.write(f'{snv+1} {cell+1} {AD.iloc[snv, cell]}\n')

        gt_file = f'{out_file_raw[0]}_gt.csv'
        pd.DataFrame(gt4, index=loci, columns=cells).to_csv(gt_file)

        DP_file = f'{out_file_raw[0]}_DP.csv'
        DP.T.to_csv(DP_file, index=True, header=True, index_label='cell_id')

        AD_file = f'{out_file_raw[0]}_AD.csv'
        AD.to_csv(AD_file, index=True, header=True, index_label='cell_id')

        RD_file = f'{out_file_raw[0]}_RD.csv'
        RD.to_csv(RD_file, index=True, header=True, index_label='cell_id')


def save_vcf(df, out_file_base):
    def conv(x):
        RD, AD, GT = x.split(':')
        if GT == '3':
            return '.:.:.'
        else:
            return f'{GT}:{int(RD)+int(AD)}:{AD}'

    cont = ''
    for i in df['CHR'].unique():
        cont += f'##contig=<ID={i}>\n'

    df.iloc[:,7:] = df.iloc[:,7:].applymap(conv)
    body = ''
    for _, row in df.iterrows():
        body += VCF_BODY_ROW.format(
            chr=row['CHR'],
            pos=row['POS'],
            ref=row['REF'],
            alt=row['ALT'],
            data='\t'.join(row[7:]),
        )
    with open(f'{out_file_base}.vcf', 'w') as f:
        f.write(VCF_HEADER.format(contigs=cont, cells='\t'.join(df.columns[7:])))
        f.write(body)


def filter_variants(df, gt, VAF, args):
    ms1 = (gt == 0) & (VAF > args.max_ref_VAF)
    ms2 = (gt == 1) & (VAF < args.min_het_VAF)
    ms3 = (gt == 2) & (VAF < args.min_hom_VAF)
    ms = ms1 | ms2 | ms3

    gt[ms] = 3
    df.iloc[:,7:] = df.iloc[:,7:] \
        .where(~ms, df.iloc[:,7:].applymap(lambda x: x[:-1] + '3'))

    keep_var1 = np.mean(gt == 3, axis=1) < (1 - args.minVarGeno)
    if args.minMutated < 1:
        keep_var2 = np.mean((gt == 1) | (gt == 2), axis=1) >= args.minMutated
    elif args.minMutated == 1:
        raise TypeError('--minMutated cannot be equactly 1. Values <1 are ' \
            'interpreted as cell fraction, values >1 as absolute cell number.')
    else:
        keep_var2 = np.sum((gt == 1) | (gt == 2), axis=1) >= args.minMutated

    keep = keep_var1 & keep_var2
    keep[get_WL_ids(df)] = True

    return df[keep].reset_index(drop=True), gt[keep]


def filter_variants_consecutive(df, gt, proximity):
    keep = np.ones(df.shape[0], dtype=bool)

    chrom = df['CHR'].values
    pos = df['POS'].values

    loc = 0
    while loc < len(chrom):
        found = 0
        fa = np.argwhere((chrom == chrom[loc]) & (np.arange(len(chrom)) > loc)) \
            .flatten()[::-1]
        for jj in fa:
            for ii in range(len(proximity)):
                if (pos[loc] + proximity[ii] > pos[jj]) & (jj - loc > ii):
                    found = 1
            if found == 1:
                keep[np.arange(loc, jj + 1)] = False
                loc = jj + 1
                break
        if found == 0:
            loc += 1

    keep[get_WL_ids(df)] = True
    return df[keep].reset_index(drop=True), gt[keep]


def filter_variants_WL(df, gt, args):
    keep = np.ones(df.shape[0], dtype=bool)
    uniq_vars = df[df.duplicated(subset=['CHR', 'POS'])][['CHR', 'POS']]

    for _, (chrom, pos) in uniq_vars.iterrows():
        var_ids = df[(df['CHR'] == chrom) & (df['POS'] == pos)].index.values
        freq = ((gt[var_ids] == 1) | (gt[var_ids] == 2)).mean(axis=1)
        bases = df.loc[var_ids, ['REF', 'ALT']]

        # Both reported ALT have a 0 frequency, display only first one
        if freq.max() == 0:
            keep[var_ids[1:]] = False
        # All or all but one low frequent, only kept cause of whitelist
        elif freq.max() < args.minMutated or sum(freq >= args.minMutated) == 1:
            keep[var_ids[freq != freq.max()]] = False
        # All but one variant have a '*' as REF or ALT
        elif ((bases == '*').sum(axis=1) == 0).sum() == 1:
            keep[var_ids[((bases == '*').sum(axis=1) > 0)]] = False
        else:
            # Keep everything but variants with '*'
            keep[var_ids[((bases == '*').sum(axis=1) > 0)]] = False
            # import pdb; pdb.set_trace()
    return df[keep].reset_index(drop=True), gt[keep]


def postprocess_mosaic(df, gt):
    rmv_all = []

    for ampl_sgl, ampl_data_raw in df.groupby('NAME'):
        # Filter variants with '*' as REF or ALT
        bad_q = (ampl_data_raw['REF'] == '*') | (ampl_data_raw['ALT'] == '*')
        if bad_q.any():
            print('Removed bad quality variants:\n' \
                f'{ampl_data_raw.loc[bad_q].iloc[:,:5]}\n')
            rmv_all.extend(bad_q.index.values)
            ampl_data_raw = ampl_data_raw[~bad_q]

        if ampl_data_raw.shape[0] < 2:
            continue

        ampl_data = ampl_data_raw.sort_values('POS')
        # Filter multiple variants at the same position if they are shady
        dupl = ampl_data['POS'].duplicated(keep=False)
        if dupl.any():
            df_dupl = ampl_data.loc[dupl]
            if dupl.sum() > 2:
                print('Removed variants with >2 ALT at the same position:\n' \
                    f'{ampl_data.loc[dupl].iloc[:,:5]}\n')
                rmv_all.extend(df_dupl.index.values)
                continue
            # Remove if REF0 = ALT1 and REF1 = ALT0 in SAME POSITION
            REF0, ALT0 = df_dupl.iloc[0, 2:4]
            REF1, ALT1 = df_dupl.iloc[1, 2:4]
            if (REF0 == ALT1) & (REF1 == ALT0):
                print(f'Removed shady alignment (unclear indel):\n'
                    f'{df_dupl.iloc[:,:5]}\n')
                rmv_all.extend(df_dupl.index.values)
                ampl_data = ampl_data_raw[~dupl]
            elif len(REF0) != len(REF1) or len(ALT0) != len(ALT1):
                print(f'Removed shady variants (different length indels):\n'
                    f'{df_dupl.iloc[:,:5]}\n')
                rmv_all.extend(df_dupl.index.values)
                ampl_data = ampl_data_raw[~dupl]

        # Remove variants closer than 5 bp together;
        #    Assumes the data to be sorted by position!
        rmv_close = np.argwhere(np.diff(ampl_data['POS']) < 5).flatten()
        if rmv_close.size > 0:
            rmv_close = np.unique(np.concatenate([rmv_close, rmv_close + 1]))
            print('Removed shady alignment (too close):' \
                    f'\n{ampl_data.iloc[rmv_close,:5]}\n')
            rmv_all.extend(ampl_data.iloc[rmv_close].index.values)
            keep_idx = [i for i in range(ampl_data.shape[0]) if i not in rmv_close]
            ampl_data = ampl_data.iloc[keep_idx]

    print(f'Removed {len(rmv_all)} additional variants')
    keep = np.where(np.isin(np.arange(df.shape[0]), rmv_all), False, True)
    keep[get_WL_ids(df)] = True

    return df[keep].reset_index(drop=True), gt[np.argwhere(keep).flatten()]


def concat_str_arrays(arrays, sep='_'):
    out_arr = arrays[0].astype(str)
    for arr in arrays[1:]:
        out_arr = np.char.add(out_arr, sep)
        out_arr = np.char.add(out_arr, arr.astype(str))
    return out_arr


def merge_gt(gt_in):
    # 00=0, 11=1, 22=2, 33=3; 03=0, 13=1, 23=2; 01=1, 02=1, 12=1
    if gt_in[0] == gt_in[1]:
        return gt_in[0]
    elif gt_in[0] == 3:
        return gt_in[1]
    elif gt_in[1] == 3:
        return gt_in[0]
    else:
        return 1


def multiplex_looms(args):
    no_samples = len(args.input)
    assert no_samples == len(args.ratio), \
        'No. input files has to be the same as no. ratios. ' \
            f'({no_samples} != {len(args.ratio)})'
    assert np.sum(args.ratio) <= 1, 'Ratios cannot sum up to >1.'

    no_cells = np.zeros(no_samples, dtype=int)
    for i, in_file in enumerate(args.input):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_in_file = os.path.join(temp_dir, os.path.basename(in_file))
            shutil.copy2(in_file, temp_in_file)
            with loompy.connect(temp_in_file) as ds:
                no_cells[i] = ds.col_attrs['barcode'].size

    # Subsample cells
    samples_total = (no_cells.min() * np.array(args.ratio)).astype(int)
    dbt_total = int(samples_total.sum() * args.doublets)
    samples_size = samples_total + (dbt_total * np.array(args.ratio)).astype(int)

    samples = {}
    for i, size in enumerate(samples_size):
        cells_multiplex = np.random.choice(no_cells[i], size=size, replace=False)
        samples[i] = {'idx': np.sort(cells_multiplex)}

    for i, in_file in enumerate(args.input):
        print(f'Reading: {in_file}')
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_in_file = os.path.join(temp_dir, os.path.basename(in_file))
            shutil.copy2(in_file, temp_in_file)

            # Open loom file and read data
            with loompy.connect(temp_in_file) as ds:
                index = concat_str_arrays(
                    [ds.ra['CHROM'], ds.ra['POS'], ds.ra['REF'], ds.ra['ALT']])
                cols = samples[i]['idx'] + ((i + 1) / 10)

                df_new = pd.DataFrame(ds[:, samples[i]['idx']],
                    index=index, columns=cols)
                ampl_new = pd.Series(ds.ra['amplicon'], index=index, name=0)
                DP_new = pd.DataFrame(ds.layers['DP'][:,samples[i]['idx']],
                    index=index, columns=cols)
                GQ_new = pd.DataFrame(ds.layers['GQ'][:,samples[i]['idx']],
                    index=index, columns=cols)
                AD_new = pd.DataFrame(ds.layers['AD'][:,samples[i]['idx']],
                    index=index, columns=cols)
                RO_new = pd.DataFrame(ds.layers['RO'][:,samples[i]['idx']],
                    index=index, columns=cols)

                barcodes = np.char.add(
                    ds.col_attrs['barcode'][samples[i]['idx']], f'.pat{i:.0f}')
                barcode_length = int(barcodes.dtype.str[2:])
                samples[i]['name'] = barcodes.astype(f'<U{2*barcode_length+1}')

        # First sample, nothing to merge
        if i == 0:
            df = df_new
            ampl = ampl_new
            DP = DP_new
            GQ = GQ_new
            AD = AD_new
            RO = RO_new
            continue

        # Merge amplicons (keep all)
        ampl = ampl.combine_first(ampl_new)
        df = df.merge(df_new, left_index=True, right_index=True)
        DP = DP.merge(DP_new, left_index=True, right_index=True)
        GQ = GQ.merge(GQ_new, left_index=True, right_index=True)
        AD = AD.merge(AD_new, left_index=True, right_index=True)
        RO = RO.merge(RO_new, left_index=True, right_index=True)

    del df_new
    del ampl_new
    del DP_new
    del GQ_new
    del AD_new
    del RO_new
    gc.collect()

    # Add doublets (if doublets arg specified)
    print('Generating doublets')
    samples['dbt'] = {'name': []}
    s_probs = np.array(args.ratio) / sum(args.ratio)

    drop_cells = []
    dbt_gt = {}
    dbt_DP = {}
    dbt_GQ = {}
    dbt_AD = {}
    dbt_RO = {}
    for i in range(dbt_total):
        s1, s2 = np.random.choice(no_samples, size=2, replace=False, p=s_probs)

        c1_idx = np.random.choice(samples[s1]['idx'].size)
        c2_idx = np.random.choice(samples[s2]['idx'].size)

        c1 = samples[s1]['idx'][c1_idx]  + ((s1 + 1) / 10)
        c2 = samples[s2]['idx'][c2_idx]  + ((s2 + 1) / 10)

        new_id = f'{c1}+{c2}'
        new_name = f'{samples[s1]["name"][c1_idx]}+{samples[s2]["name"][c2_idx]}'

        dbt_gt[new_id] = np.apply_along_axis(merge_gt, axis=1, arr=df[[c1, c2]])
        dbt_GQ[new_id] = GQ[[c1,c2]].mean(axis=1)
        dbt_DP[new_id] = (DP[[c1,c2]].mean(axis=1).round()).astype(int)
        dbt_AD[new_id] = (AD[[c1,c2]].mean(axis=1).round()).astype(int)
        dbt_RO[new_id] = dbt_DP[new_id] - dbt_AD[new_id]

        samples[s1]['idx'] = np.delete(samples[s1]['idx'], c1_idx)
        samples[s1]['name'] = np.delete(samples[s1]['name'], c1_idx)
        samples[s2]['idx'] = np.delete(samples[s2]['idx'], c2_idx)
        samples[s2]['name'] = np.delete(samples[s2]['name'], c2_idx)
        samples['dbt']['name'].append(new_name)

        drop_cells.extend([c1, c2])

    df = df.join(pd.DataFrame(dbt_gt, index=df.index))
    DP = DP.join(pd.DataFrame(dbt_DP, index=df.index))
    GQ = GQ.join(pd.DataFrame(dbt_GQ, index=df.index))
    AD = AD.join(pd.DataFrame(dbt_AD, index=df.index))
    RO = RO.join(pd.DataFrame(dbt_RO, index=df.index))

    df.drop(drop_cells, axis=1, inplace=True)
    DP.drop(drop_cells, axis=1, inplace=True)
    GQ.drop(drop_cells, axis=1, inplace=True)
    AD.drop(drop_cells, axis=1, inplace=True)
    RO.drop(drop_cells, axis=1, inplace=True)

    gt = df.fillna(3).values
    index = df.index.values
    cells = df.columns.values

    del df
    gc.collect()

    print('Filtering data')
    mut = (gt == 1) | (gt == 2)
    if args.minMutated < 1:
        mut_var = np.mean(mut, axis=1) >= args.minMutated
    else:
        mut_var = np.sum(mut, axis=1) >= args.minMutated

    gt = gt[mut_var]
    DP = DP[mut_var].values
    GQ = GQ[mut_var].values
    AD = AD[mut_var].values
    RO = RO[mut_var].values
    index = index[mut_var]
    ampl = ampl.loc[index].values
    # Filters 1-3: done second in mosaic
    # First filter: Remove genotype in cell with quality < X
    gt[GQ < args.minGQ] = 3
    del GQ
    # Second filter: Remove genotype in cell with read depth < X
    gt[DP < args.minDP] = 3
    # Third filter: Remove genotype in cell with alternate allele freq < X
    VAF = np.where(DP > 0, AD / DP, 0)
    gt[((gt == 1) | (gt == 2)) & (VAF < args.minVAF)] = 3

    del DP
    gc.collect()
    # Fourth filter: Remove variants genotyped in < X % of cells; done last/fourth in mosaic
    keep_var1 = np.mean(gt == 3, axis=1) < (1 - args.minVarGeno)
    # Fifth filter: Remove cells with < X %of genotypes present; done third in mosaic
    keep_cells = np.mean(gt[keep_var1] == 3, axis=0) < (1 - args.minCellGeno)

    # keep_cells = np.full(cells.size, True)
    cell_names = []
    for i in samples:
        cell_names.extend(samples[i]['name'])

    cells = cells[keep_cells]
    cell_names = np.array(cell_names)[keep_cells]
    # Filter second time for variants present in >X percent of data
    gt = gt[:,keep_cells]
    AD = AD[:,keep_cells]
    RO = RO[:,keep_cells]
    VAF = VAF[:,keep_cells]

    if args.minMutated < 1:
        keep_var2 = np.mean((gt == 1) | (gt == 2), axis=1) >= args.minMutated
    else:
        keep_var2 = np.sum((gt == 1) | (gt == 2), axis=1) >= args.minMutated

    # Get only variants passing both filters
    keep_var = keep_var1 & keep_var2

    gt = gt[keep_var]
    index = index[keep_var].astype(str)
    ampl = ampl[keep_var]
    RO = RO[keep_var]
    AD = AD[keep_var]
    VAF = VAF[keep_var]

    meta_info = np.stack(np.char.split(index, '_'))
    variants_info = {
        'CHR': meta_info[:,0],
        'POS': meta_info[:,1].astype(int),
        'REF': meta_info[:,2],
        'ALT': meta_info[:,3],
        'REGION': ampl,
        'NAME': ampl,
        'FREQ': np.zeros(np.sum(keep_var))
    }

    for row, var_data in enumerate(gt):
        for col, cell in enumerate(cells):
            name = cell_names[col]
            cell_str = f'{RO[row, col]}:{AD[row, col]}:{var_data[col]}'
            try:
                variants_info[name].append(cell_str)
            except KeyError:
                variants_info[name] = [cell_str]

    return pd.DataFrame(variants_info), gt, VAF



def process_loom(args, rel_cells=np.array([])):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_in_file = os.path.join(temp_dir, os.path.basename(args.input[0]))
        shutil.copy2(args.input[0], temp_in_file)
        with loompy.connect(temp_in_file) as ds:
            gt = ds[:,:]
            if rel_cells.size > 0:
                cell_names = ds.col_attrs['barcode']
                cells = np.argwhere(np.isin(cell_names, rel_cells)).flatten()
            else:
                # Subsample cells
                _, cell_no = gt.shape
                if np.isfinite(args.cell_no):
                    if int(args.cell_no) == args.cell_no:
                        subsample_size = min(int(args.cell_no), cell_no)
                    else:
                        subsample_size = min(int(args.cell_no * cell_no), cell_no)
                    print(f'Subsampling {subsample_size} cells')
                    cells = np.random.choice(cell_no, size=subsample_size, replace=False)
                    cells = np.sort(cells)
                else:
                    cells = np.arange(cell_no)

            gt = gt[:, cells]

            # Sixth filter: Remove variants mutated in < X % of cells; done first in mosaic
            mut = (gt == 1) | (gt == 2)

            if args.minMutated < 1:
                mut_var = np.mean(mut, axis=1) >= args.minMutated
            else:
                mut_var = np.sum(mut, axis=1) >= args.minMutated

            # ids = np.char.add(np.char.add(ds.ra['CHROM'], '_'), ds.ra['POS'].astype(str))
            var_id_all = concat_str_arrays([ds.ra['CHROM'], ds.ra['POS']])
            wl_var_id = np.argwhere(np.isin(var_id_all, WHITELIST)).flatten()

            mut_var[wl_var_id] = True

            gt = gt[mut_var]
            var_id = var_id_all[mut_var]

            DP = ds.layers['DP'][mut_var,:][:,cells]
            GQ = ds.layers['GQ'][mut_var,:][:,cells]
            AD = ds.layers['AD'][mut_var,:][:,cells]
            RO = ds.layers['RO'][mut_var,:][:,cells]

            # Filters 1-3: done second in mosaic
            # First filter: Remove genotype in cell with quality < X
            gt[GQ < args.minGQ] = 3
            # Second filter: Remove genotype in cell with read depth < X
            gt[DP < args.minDP] = 3
            # Third filter: Remove genotype in cell with alternate allele freq < X
            VAF = np.where(DP > 0, AD / DP, 0)
            gt[((gt == 1) | (gt == 2)) & (VAF < args.minVAF)] = 3

            del DP
            gc.collect()
             # Fourth filter: Remove variants genotyped in < X % of cells; done last/fourth in mosaic
            keep_var1 = np.mean(gt == 3, axis=1) < (1 - args.minVarGeno)
            # Fifth filter: Remove cells with < X %of genotypes present; done third in mosaic
            keep_cells = np.mean(gt[keep_var1] == 3, axis=0) < (1 - args.minCellGeno)

            # Filter second time for variants present in >X percent of data
            gt = gt[:,keep_cells]
            keep_var2 = np.mean((gt == 1) | (gt == 2), axis=1) > 0.01
            # Get only variants passing both filters
            keep_var = keep_var1 & keep_var2
            keep_var_wl = np.argwhere(np.isin(var_id, WHITELIST)).flatten()

            keep_var[keep_var_wl] = True
            gt = gt[keep_var]

            RO = RO[keep_var][:,keep_cells]
            AD = AD[keep_var][:,keep_cells]
            VAF = VAF[keep_var][:,keep_cells]

            ampl = ds.ra['amplicon'][mut_var][keep_var]
            
            variants_info = {
                'CHR': ds.ra['CHROM'][mut_var][keep_var],
                'POS': ds.ra['POS'][mut_var][keep_var],
                'REF': ds.ra['REF'][mut_var][keep_var],
                'ALT': ds.ra['ALT'][mut_var][keep_var],
                'REGION': ampl,
                'NAME': ampl,
                'FREQ': np.zeros(np.sum(keep_var))
            }
            cells = ds.col_attrs['barcode'][cells][keep_cells]

    for cell in cells:
        variants_info[cell] = []

    for row, var_data in enumerate(gt):
        for col, cell in enumerate(cells):
            cell_str = f'{RO[row, col]}:{AD[row, col]}:{var_data[col]}'
            try:
                variants_info[cell].append(cell_str)
            except KeyError:
                variants_info[cell] = [cell_str]

    return pd.DataFrame(variants_info), gt, VAF


def update_whitelist(wl_file):
    df = pd.read_csv(wl_file, sep='\t')

    if df.size == 0:
        if df.shape[1] > 1:
            for snv in df.columns:
                WHITELIST.append('_'.join(snv.split('_')[:2]))
        else:
            for snv in df.columns.values[0].split(';'):
                WHITELIST.append('_'.join(snv.split('_')[:2]))
    # Just 1 column: either wrong separator or direct whitelist entries
    elif df.shape[1] == 1:
        if re.match('[\dXY]{1,2}_\d+', df.columns[0]):
            with open(wl_file, 'r') as f:
                ids = f.read().strip().split('\n')
        else:
            df = pd.read_csv(wl_file, sep=',')
            ids = df['CHR'].astype(str) + '_' + df['POS'].astype(str)
        WHITELIST.extend(ids)
    elif 'chr' in df.columns and 'pos' in df.columns:
        for _, (chrom, pos) in df.iterrows():
            WHITELIST.append(f'{chrom.replace("chr", "")}_{pos}')
    else:
        for _, (sample, snvs) in df.iterrows():
            if '+' in sample:
                continue
            for snv in snvs.split(';'):
                WHITELIST.append('_'.join(snv.split('_')[:2]))


def parse_args():
    parser = argparse.ArgumentParser(prog='mosaic preprocessing',
        usage='Single loom: python mosaic_preprocessing.py -i <DATA> [-args];' \
            'Multiplexing: python mosaic_preprocessing.py -i <LOOM_1> <LOOM_2> ' \
            '<LOOM_X> -R <RATIO_1> <RATIO_2> <RATIO_X> -d <DOUBLET_RATE>',
        description='*** Filter a single loom file equal to MissionBio mosaic, ' \
        'or multiplex n loom files synthetically and filter afterwards. ***')
    parser.add_argument('-i', '--input', nargs='+', type=str,
        help='Input loom file. If more given, multiplex snythetically.')
    parser.add_argument('-a', '--assignment', type=str, default='',
        help='Assignment of cells to clusters. Create output for each cluster.')
    parser.add_argument('-pa', '--panel', type=str, default='',
        help='Path to (panel) Designer folder with XX-submitted.bed and ' \
        'XX-amplicon.bed file. Required for defining genes as regions.')
    parser.add_argument('-wl', '--whitelist', type=str, default='',
        help='Whitelist containing SNVs for plotting. Default = None.')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('-fo', '--full_output', action='store_true',
        help='Generate more output files (AD, RD, DP, full GT).')

    multiplex = parser.add_argument_group('multiplex simulations')
    multiplex.add_argument('-r', '--ratio', nargs='+', type=float,
        help='Multiplex ratio (#arguments == #inputs).')
    multiplex.add_argument('-d', '--doublets', type=float, default=0,
        help='Fraction of doublet. Default = 0.')

    downsample = parser.add_argument_group('downsample simulations')
    downsample.add_argument('-c', '--cell_no', type=float, default=np.inf,
        help='Downsample cells (either total number or fraction).' \
            'Default = None')

    mosaic = parser.add_argument_group('mosaic')
    mosaic.add_argument('-mStd', '--minStd', type=float, default=0.1,
        help='Minimum StandardDeviation at a loci to consider a variant.' \
            'If the Std is lower, the loci is considered a germline/artefact')
    mosaic.add_argument('-mGQ', '--minGQ', type=int, default=30,
        help='Minimum GQ to consider a variant (per cell).')
    mosaic.add_argument('-mDP', '--minDP', type=int, default=10,
        help='Minimum depth to consider a variant (per cell) .')
    mosaic.add_argument('-mAF', '--minVAF', type=float, default=0.2,
        help='Minimum alternate VAF to consider a variant (per cell).')
    mosaic.add_argument('-mVG', '--minVarGeno', type=float, default=0.5,
        help='Minimum fractions of genotyped loci to consider a variant.')
    mosaic.add_argument('-mCG', '--minCellGeno', type=float, default=0.5,
        help='Minimum fractions of genotyped cells to consider a cell.')
    mosaic.add_argument('-mMut', '--minMutated', type=float, default=0.01,
        help='Minimum fractions or cell number of het/hom cells to consider a variant.')
    mosaic.add_argument('-refVAF', '--max_ref_VAF', type=float, default=0.05,
        help='Maximum VAF to consider a wildtype variant (per cell).')
    mosaic.add_argument('-homVAF', '--min_hom_VAF', type=float, default=0.95,
        help='Minimum VAF to consider a homozygous variant (per cell).')
    mosaic.add_argument('-hetVAF', '--min_het_VAF', type=float, default=0.35,
        help='Minimum VAF to consider a heterozygous variant (per cell).')
    mosaic.add_argument('-p', '--proximity', nargs='+', type=int,
        default=[25, 50, 100, 200], help='If "i + 1" variants are within ' \
            '"proximity[i]", then the variants are removed.')
    mosaic.add_argument('-fs', '--filter_same', action='store_true',
        help='Filter all but one variant on the same amplicon ' \
            'if their genotypes are >95%% equal. Default = False')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.whitelist:
        update_whitelist(args.whitelist)
    main(args)
