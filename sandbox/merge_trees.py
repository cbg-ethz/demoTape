#!/usr/bin/env python3

import argparse
import os
import re


class Tree:
    STYLE = '[color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5]'
    GV = 'digraph G{{\nnode {style};\n{mut_edges}\n{mut_nodes}\n{cell_edges}\n{cell_nodes}\n}}'


    def __init__(self):
        self.nodes = {}
        self.edges = []


    def read(self, tree_str):
        # Init nodes
        for line in tree_str.split('\n'):
            if line.startswith(('digraph', 'node', '}')) or '->' in line:
                continue
            break_idx = line.index('[')
            name = line[:break_idx].strip()
            style = line[break_idx:].rstrip(';')
            self.nodes[name] = Node(name, style)

        # Init edges
        for line in tree_str.split('\n'):
            if line.startswith(('digraph', 'node', '}')) or not '->' in line:
                continue

            break_idx = line.index('[')
            start_id, end_id = line[:break_idx].split('->')
            start = self.nodes[start_id.strip()]
            end = self.nodes[end_id.strip()]
            end.not_root()
            style = line[break_idx:].rstrip(';')
            self.edges.append(Edge(start, end, style))

        self.get_root() # Sanity check for 1 root node


    def combine(self, trees):
        new_root_name = '0'
        new_root = Node(new_root_name, '[label=<MERGE ROOT<br/>>]')
        self.nodes[new_root_name] = new_root
        for tree in trees:
            self.nodes.update(tree.nodes)
            self.edges.extend(tree.edges)

            old_root = tree.get_root()
            old_root.not_root() # anymore
            self.edges.append(
                Edge(new_root, old_root, '[color=dimgray penwidth=4 weight=2]'))


    def no_nodes(self, n_type='all'):
        if n_type == 'all':
            return len(self.nodes)
        else:
            return len([i for i in self.nodes.values() if i.type == n_type])


    def increase_node_id(self, mut_add, cell_add):
        new_nodes = {}
        for node_id, node in self.nodes.copy().items():
            old_id = node_id
            if node.type == 'mut':
                new_id = str(int(node_id) + mut_add)
            else:
                new_id = str(int(node_id) + cell_add)
            node.set_name(new_id)
            new_nodes[new_id] = self.nodes[old_id]
        self.nodes = new_nodes


    def get_root(self):
        roots = [i for i in self.nodes.values() if i.is_root]
        assert len(roots) == 1, '!= 1 root node detected'
        return roots[0]


    def get_nodes(self, n_type='all'):
        if n_type == 'all':
            return [i for i in self.nodes.values()]
        else:
            return [i for i in self.nodes.values() if i.type == n_type]


    def get_edges(self, e_type='all'):
        if e_type == 'all':
            return self.edges
        else:
            return [i for i in self.edges if i.type == e_type]


    def __str__(self):
        mut_nodes = '\n'.join([str(i) for i in self.get_nodes('mut')])
        mut_edges = '\n'.join([str(i) for i in self.get_edges('mut')])
        cell_nodes = '\n'.join([str(i) for i in self.get_nodes('cell')])
        cell_edges = '\n'.join([str(i) for i in self.get_edges('cell')])
        return Tree.GV.format(
            style=Tree.STYLE,
            mut_nodes=mut_nodes,
            mut_edges=mut_edges,
            cell_nodes=cell_nodes,
            cell_edges=cell_edges)


class Node:
    def __init__(self, name, style, root=True):
        self.name = name
        self.style = style

        if 'style' in style:
            self.type = 'cell'
            self.is_root = False
        else:
            self.type = 'mut'
            self.is_root = root
            events = style[8:-7].split('<br/>')
            SNVs = [i for i in events if not i.startswith('<B>') if i]
            # CNVs = [i for i in events if i.startswith('<B>')]
            for SNV in SNVs:
                if SNV == 'MERGE ROOT':
                    continue
                if not re.search('.*\(chr[0-9XY]+_\d+\)', SNV):
                    print(f'!WARNING: mutation event "{SNV}" does not contain ' \
                        ' position (format: chr<CHR>_<POS>)')

    def __str__(self):
        return f'{self.name} {self.style};'


    def set_name(self, name):
        self.name = name


    def not_root(self):
        self.is_root = False


class Edge:
    def __init__(self, start, end, style):
        self.start = start
        self.end = end
        self.style = style.rstrip(';')
        if 'dir=none' in style:
            self.type = 'cell'
        else:
            self.type = 'mut'


    def __str__(self):
        return f'{self.start.name} -> {self.end.name} {self.style};'


def merge_trees(trees):
    # Get node type specific numbers
    mut_nodes = [i.no_nodes('mut') for i in trees]
    cell_nodes = [i.no_nodes('cell') for i in trees]

    # Adjust node indices: mutation nodes 1, ..., n; cell nodes n + 1, ..., m
    for i, tree in enumerate(trees):
        # +1 for additional root node
        tree.increase_node_id(
            mut_add=sum(mut_nodes[:i]) + 1,
            cell_add=sum(mut_nodes) + 1 + sum(cell_nodes[:i]) - mut_nodes[i])

    out_tree = Tree()
    out_tree.combine(trees)
    return out_tree



def main(args):
    trees = []
    for in_file in args.input:
        with open(in_file, 'r') as f:
            tree_str = f.read().strip()
        new_tree = Tree()
        new_tree.read(tree_str)
        trees.append(new_tree)

    out_tree = merge_trees(trees)

    if not args.output:
        args.output = args.input[0] + '.merged.gv'

    with open(args.output, 'w') as f:
        f.write(str(out_tree))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, nargs='+',
        help='Input .gz tree files.')
    parser.add_argument('-o', '--output', type=str, help='Output tree file.')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)