#!/usr/bin/env python3

import argparse

import matplotlib.pyplot as plt
import missionbio.mosaic as ms
import numpy as np


MARGINS = {
    'left': 0.15,
    'right': 0.9,
    'top': 0.9,
    'bottom': 0.1,
    'wspace': 0.25,
    'hspace': 0.75,
}


def main(args):
    sample = ms.load(args.input, raw=False, apply_filter=True, single=True)
    filtered_variants1 = sample.dna.filter_variants()
    print(f'# Vars before filtering: {sample.dna.shape[1]}')
    sample.dna = sample.dna[sample.dna.barcodes(), filtered_variants1]
    print(f'# Vars  after 1. filtering: {sample.dna.shape[1]}')
    filtered_variants2 = sample.dna.filter_variants_consecutive()
    sample.dna = sample.dna[sample.dna.barcodes(), filtered_variants2]
    print(f'# Vars  after 2. filtering: {sample.dna.shape[1]}')

    annotation = sample.dna.get_annotations()
    for col, content in annotation.items():
        sample.dna.add_col_attr(col, content.values)

    clusters = sample.dna.group_by_genotype(features=filtered_variants2)
    sample.dna.add_row_attr("clustering-1", sample.dna.get_labels())
    print(np.unique(sample.dna.get_labels(), return_counts=True))

    sample.dna.find_clones()
    sample.dna.add_row_attr("clustering-2", sample.dna.get_labels())
    print(np.unique(sample.dna.get_labels(), return_counts=True))

    rel_variants = sample.dna.find_relevant_variants()
    print(f'# Relevant Vars: {len(rel_variants)}')

    scatter = sample.dna.scatterplot(attribute='umap', colorby='label')
    if not args.output:
        scatter.show()
    else:
        scatter.write_image(args.output + '.scatter.png')

    hm = sample.dna.heatmap(attribute='AF')
    if not args.output:
        hm.show()
    else:
        hm.write_image(args.output + '.heatmap.png')

    import pdb; pdb.set_trace()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input h5 file.')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)