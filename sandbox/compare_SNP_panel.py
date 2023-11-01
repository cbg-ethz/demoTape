#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd

DEFAULT_PANEL = '/home/hw/Desktop/mincedTapestri/resources/Designer/3848-amplicon.bed'


def map_snp_to_panel(snp_file, panel_file):
    panel = pd.read_csv(panel_file, sep='\t', header=None,
        names=['chr', 'start', 'stop', 'amplicon'], skiprows=[0], index_col=0)

    try:
        # SNP array profiles
        snps = pd.read_csv(snp_file, sep='\t')
        snp_type = 'SNP array'
    except pd.errors.ParserError:
        # scRNAseq profiles
        snps = pd.read_csv(snp_file, sep='\t', comment='#', header=None,
            usecols=[0, 1, 3, 4, 6],
            names=['#CHROM', 'POS', 'REF', 'ALT', 'FILTER'])
        snps = snps[(snps['FILTER'] == '.') | (snps['FILTER'] == 'PASS')]
        snp_type = 'scRNAseq SNPs'

    if snps.iloc[0, 0].startswith('chr'):
        snps['#CHROM'] = snps['#CHROM'].str.replace('chr', '')
    snps['name'] = snps['#CHROM'] + '_' + snps['POS'].astype(str)

    snp_map = {}
    print('Mapping SNP file to panel')
    for chrom, snp_chr_df in snps.groupby('#CHROM'):
        try:
            p_chr_df = panel.loc[f'chr{chrom}']
        except KeyError:
            continue

        overlap = snp_chr_df[(snp_chr_df['POS'] > p_chr_df['start'].min()) \
            & (snp_chr_df['POS'] < p_chr_df['stop'].max())]

        for i, snp in overlap.iterrows():
            pos = snp['POS']
            match = (pos >= p_chr_df['start']) & (pos <= p_chr_df['stop'])
            if np.any(match):
                ampl = p_chr_df.iloc[match.to_numpy().nonzero()[0][0], 2]
                snp_map[snp['name']] = (ampl, snp['REF'], snp['ALT'])

    print(f'\nOverlap {snp_type} & panel: {len(snp_map)}')
    return snp_map


def map_variants_to_SNPs(snp_map, variants_file, verbose=False):
    variants = pd.read_csv(variants_file, index_col=[0,1])
    print('SNP panel matches:')

    map_data = []
    for snp, ampl in snp_map.items():
        chrom, pos = snp.split('_')
        try:
            var = variants.loc[chrom, int(pos)]
        except (TypeError, KeyError):
            pass
        else:
            gt = np.array(var[5:].apply(lambda x: int(x.split(':')[-1])).values)
            wt = np.mean(gt == 0)
            het = np.mean(gt == 1)
            hom = np.mean(gt == 2)
            mis = np.mean(gt == 3)
            if ampl[1] != var["REF"] or ampl[2] != var["ALT"]:
                print(
                    f'\tDifferent REF|ALT in SNP panel: {ampl[1]} -> {ampl[2]}' \
                    f'\n\tscDNAseq: chr{snp} ({var["REF"]} -> {var["ALT"]}) - ' \
                    f'wt: {wt:.2f}, het: {het:.2f}, hom: {hom:.2f}, missing: {mis:.2f}')
            if verbose:
                print(f'\tchr{snp: <14} ({var["REF"]} -> {var["ALT"]}) - ' \
                    f'wt: {wt:.2f}, het: {het:.2f}, hom: {hom:.2f}, missing: {mis:.2f}')
            map_data.append([chrom, pos, ampl[1], ampl[2], wt, het, hom, mis])

    df = pd.DataFrame(map_data,
        columns=['CHR', 'POS', 'REF', 'ALT', 'WT', 'HET', 'HOM', 'MISSING'])

    def _to_sort_idx(chrom, pos):
        if chrom.isdigit():
            return float(f'{chrom}.{pos}')
        else:
            return float(f'23.{pos}')
    df['idx'] = df.apply(lambda x: _to_sort_idx(x['CHR'], x['POS']), axis=1)
    df.sort_values('idx', inplace=True)
    df.drop('idx', axis=1, inplace=True)
    print(f'\nMatches total: {df.shape[0]}\n')
    return df


def save_overlap(overlap, outfile):
    with open(outfile, 'w') as f:
        f.write('SNP\tAMPLICON\tREF\tALT\n')
        for snp, details in overlap.items():
            f.write(f'{snp}\t{details[0]}\t{details[1]}\t{details[2]}\n')


def main(args):
    snp_map = map_snp_to_panel(args.SNPs, args.panel)
    # save_overlap(snp_map, args.outfile)
    var_map = map_variants_to_SNPs(snp_map, args.variants, args.outfile == None)
    if args.outfile:
        var_map.to_csv(args.outfile, index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--SNPs', type=str, help='SNP .vcf array file.')
    parser.add_argument('-v', '--variants', type=str, help='Variant .csv file.')
    parser.add_argument('-p', '--panel', type=str, help='Panel .bed file.')
    parser.add_argument('-o', '--outfile', type=str, help='Output file for matches.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if not args.panel:
        args.panel = DEFAULT_PANEL
    main(args)