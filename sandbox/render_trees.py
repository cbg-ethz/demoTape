#!/usr/bin/env python3

import argparse
import os
import re
import subprocess
import tempfile

from tqdm import tqdm

CHROM = {**{str(i): i for i in range(1, 23, 1)}, **{'X': 23, 'Y': 24}}
COLORS = {'cl0': '#FF7575', 'cl1': '#84B5DE', 'cl2': '#AAE3A7',
    'S3': '#e41a1c', 'S1': '#377eb8', 'S2': '#4daf4a'}

COLORS = {'cl0': '#e41a1c', 'cl1': '#377eb8', 'cl2': '#4daf4a',
    'S3': '#FF7575', 'S1': '#84B5DE', 'S2': '#AAE3A7'}


class Tree:
    STYLE = '[color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5]'
    GV = 'digraph G{{\nnode {style};\n{mut_edges}\n{mut_nodes}\n{cell_edges}\n{cell_nodes}\n}}'


    def __init__(self, tree_file):
        self.nodes = {}
        self.edges = []
        self.file = tree_file
        self._read(tree_file)


    def _read(self, tree_file):
        self.file = tree_file
        with open(tree_file, 'r') as f:
            tree_str = f.read().strip()
        # Init nodes
        for line in tree_str.split('\n'):
            if line.startswith(('digraph', 'node', '}')) or '->' in line:
                continue
            self._add_node(line)

        # Init edges
        for line in tree_str.split('\n'):
            if line.startswith(('digraph', 'node', '}')) or not '->' in line:
                continue
            self._add_edge(line)

        self.get_root() # Sanity check for 1 root node


    def _add_node(self, line):
        break_idx = line.index('[')
        name = int(line[:break_idx].strip())
        style = line[break_idx:].rstrip(';')

        if 'style' in style:
            self.nodes[name] = CellNode(name, style)
        else:
            self.nodes[name] = EventNode(name, style)


    def _add_edge(self, line):
        break_idx = line.index('[')
        start_id, end_id = line[:break_idx].split('->')
        start = self.nodes[int(start_id.strip())]
        end = self.nodes[int(end_id.strip())]
        if isinstance(end, EventNode):
            end.not_root()
        style = line[break_idx:].rstrip(';')

        if 'dir=none' in style:
            self.edges.append(CellEdge(start, end, style))
        else:
            self.edges.append(EventEdge(start, end, style))


    def _adjust_indices(self):
        mut_nodes = self.no_nodes('events')
        for i, node in enumerate(sorted(self.get_nodes('event'), key=lambda x: x.name)):
            node.set_name(i)

        for i, node in enumerate(sorted(self.get_nodes('cell'), key=lambda x: x.name)):
            node.set_name(i + mut_nodes)

        self.nodes = {i.name: i for i in self.get_nodes()}


    def remove_sparse_nodes(self, threshold=1):
        # Repeat until all tip nodes have > threshold cells assigned
        events = []
        while True:
            del_edges = []
            for edge in self.get_edges('cell'):
                if edge.end.get_percentage() <= threshold:
                    self.nodes.pop(edge.end.name)
                    del_edges.append(edge)

                    if self.is_tip(edge.start):
                        self.nodes.pop(edge.start.name)
                        upstr_edge = self.get_upstream_edge(edge.start)
                        upstr_cells = self.get_downstream_edges(upstr_edge.start, 'cell')
                        upstr_cells.end.add_cells(edge.end.get_cells())
                        upstr_cells.end.add_percentage(edge.end.get_percentage())

                        del_edges.append(upstr_edge)

                    elif edge.start.is_root:
                        start_outg = [i for i in self.get_edges('event') \
                            if i.start == edge.start]

                        if len(start_outg) == 1:
                            merge_node = start_outg[0].end

                            merge_cells = self.get_downstream_edges(merge_node, 'cell')
                            merge_cells.end.add_cells(edge.end.get_cells())
                            merge_cells.end.add_percentage(edge.end.get_percentage())
                            merge_cells.start = edge.start

                            downstr_start_outg = self.get_downstream_edges(merge_node, 'event')

                            edge.start.add_events(merge_node.get_events())

                            for redirect_edge in downstr_start_outg:
                                if redirect_edge.end not in [i.start for i in del_edges]:
                                    redirect_edge.start = edge.start

                            del_edges.append(start_outg[0])
                            self.nodes.pop(merge_node.name)
                    else:
                        upstr_start = self.get_upstream_edge(edge.start).start
                        upstr_start_outg = [i for i in self.get_downstream_edges(upstr_start, 'event')
                            if i not in del_edges]
                        downstr_start_outg = self.get_downstream_edges(edge.start, 'event')

                        # Merge to upstream if single edge
                        if len(upstr_start_outg) == 1:
                            new_start = upstr_start_outg[0].start
                            new_start.add_events(edge.start.get_events())

                            downstr_cell = self.get_downstream_edges(new_start, 'cell').end
                            downstr_cell.add_cells(edge.end.get_cells())
                            downstr_cell.add_percentage(edge.end.get_percentage())

                            for redirect_edge in self.get_downstream_edges(edge.start, 'event'):
                                if redirect_edge.end not in [i.end for i in del_edges]:
                                    redirect_edge.start = new_start

                            del_edges.append(self.get_upstream_edge(edge.start))
                            self.nodes.pop(edge.start.name)
                        # Merge to downstream if single edge
                        elif len(downstr_start_outg) == 1:
                            new_end = downstr_start_outg[0].end
                            new_end.add_events(edge.start.get_events())

                            downstr_cell = self.get_downstream_edges(new_end, 'cell').end
                            downstr_cell.add_cells(edge.end.get_cells())
                            downstr_cell.add_percentage(edge.end.get_percentage())

                            upstr_edge = self.get_upstream_edge(edge.start)
                            downstr_start_outg[0].start = upstr_edge.start

                            del_edges.append(upstr_edge)
                            self.nodes.pop(edge.start.name)


            for del_edge in set(del_edges):
                self.edges.remove(del_edge)

            self._adjust_indices()

            if not del_edges:
                break


    def set_cell_node_color(self, new_color):
        for node in self.get_nodes('cell'):
            node.set_color(new_color)
            self.get_upstream_edge(node).set_color(new_color)


    def no_nodes(self, n_type='all'):
        if n_type == 'all':
            return len(self.nodes)
        elif n_type == 'cell':
            return len([i for i in self.nodes.values() if isinstance(i, CellNode)])
        else:
            return len([i for i in self.nodes.values() if isinstance(i, EventNode)])


    def get_root(self):
        roots = [i for i in self.nodes.values() if i.is_root]
        assert len(roots) == 1, '!= 1 root node detected'
        return roots[0]


    def get_nodes(self, n_type='all'):
        if n_type == 'all':
            return [i for i in self.nodes.values()]
        elif n_type == 'cell':
            return [i for i in self.nodes.values() if isinstance(i, CellNode)]
        else:
            return [i for i in self.nodes.values() if isinstance(i, EventNode)]


    def get_SNV_ids(self, unique=True):
        ids = []
        for event_node in self.get_nodes('event'):
            ids.extend([i.get_pos_id() for i in event_node.get_events('SNV')])
        if unique:
            return sorted(list(set(ids)), key=lambda x: CHROM[x[0]])
        else:
            return sorted(ids, key=lambda x: CHROM[x[0]])


    def get_LOH_amplicons(self):
        ids = []
        for event_node in self.get_nodes('event'):
            ids.extend([i.ampl for i in event_node.get_LOHs()])
        return ids


    def get_edges(self, e_type='all'):
        if e_type == 'all':
            return self.edges
        elif e_type == 'cell':
            return [i for i in self.edges if isinstance(i, CellEdge)]
        else:
            return [i for i in self.edges if isinstance(i, EventEdge)]


    def get_upstream_edge(self, node):
        for edge in self.edges:
            if edge.end == node:
                return edge


    def get_downstream_edges(self, node, e_type='all'):
        edges = [i for i in self.get_edges(e_type) if i.start == node]

        if e_type == 'cell':
            try:
                return edges[0]
            except IndexError:
                return None
        else:
            return edges


    def is_tip(self, node):
        is_tip = True
        for edge in self.get_edges('event'):
            if edge.start == node:
                is_tip = False
        return is_tip


    def to_gv(self, out_file):
        with open(out_file, 'w') as f:
            f.write(str(self))


    def render(self, out_file='', file_type='png', overwrite=False):
        if not out_file:
            out_file = f'{os.path.splitext(self.file)[0]}.{file_type}'

        if os.path.exists(out_file) and overwrite:
            print(f'File existing: {out_file}')
            return

        tmp_tree_file = tempfile.NamedTemporaryFile(delete=False)
        self.to_gv(tmp_tree_file.name)



        cmmd = f'dot -T{file_type} {tmp_tree_file.name} -o {out_file}'
        dot = subprocess.Popen(cmmd,
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = dot.communicate()
        dot.wait()
        tmp_tree_file.close()

        if stderr or stdout:
            print(f'\nFAILED! Command: {cmmd}:\n')
            [print(i) for i in str(stdout).split('\\n')]
            [print(i) for i in str(stderr).split('\\n')]
            raise RuntimeError()


    def __str__(self):
        mut_nodes = '\n'.join([str(i) for i in self.get_nodes('event')])
        mut_edges = '\n'.join([str(i) for i in self.get_edges('event')])
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
        self.is_root = root
        self.sample_id = -1


    def __str__(self):
        return f'{self.name} {self.style};'


    def set_name(self, name):
        self.name = name


    def set_id(self, sample_id):
        self.sample_id = sample_id


    def get_sample_id(self):
        return self.sample_id


class EventNode(Node):
    def __init__(self, name, style, root=True):
        super().__init__(name, style, root)
        self.events = {'SNV': [], 'LOH': []}
        events = style[8:-7].split('<br/>')
        for event in events:
            if event == 'MERGE ROOT' or event == '':
                continue
            if event.startswith('<B>'):
                # Example: <B>LOH AMPL327913(chr8_59739333-ALT)</B>

                ampl = event.split('(')[0].split('LOH ')[1]
                chrom = event.split('(')[1].split('_')[0][3:]
                pos = {}
                for LOH_event in event[:-5].split('(')[1].split(','):
                    pos_id = int(LOH_event.split('-')[0].split('_')[1])
                    pos[pos_id] = LOH_event.split('-')[1]
                self.events['LOH'].append(LOH(ampl, chrom, pos))
            else:
                # Example: AMPL278226(chr3_31712170)
                try:
                    SNV_id = re.search('.*\((chr[0-9XY]+_\d+)\)', event).group(1)
                except (AttributeError, TypeError):
                    print(f'!WARNING: mutation event "{event}" does not contain ' \
                        ' position (format: chr<CHR>_<POS>)')
                else:
                    ampl = event.split('(')[0]
                    chrom = SNV_id.split('_')[0][3:]
                    pos = int(SNV_id.split('_')[1])
                    self.events['SNV'].append(SNV(ampl, chrom, pos))


    def not_root(self):
        self.is_root = False


    def get_events(self, e_type='all'):
        if e_type == 'SNV':
            return self.events['SNV']
        elif e_type == 'LOH':
            return self.events['LOH']
        else:
            return self.events['SNV'] + self.events['LOH']


    def remove_events(self, del_events):
        for e_type in self.events:
            self.events[e_type] = [i for i in self.events[e_type] \
                if i not in del_events]


    def add_events(self, new_events):
        for event in new_events:
            if isinstance(event, SNV):
                self.events['SNV'].append(event)
            elif isinstance(event, LOH):
                self.events['LOH'].append(event)
            else:
                raise TypeError(f'Unknown event: {type(event)}')


    def __str__(self, sep='<br/>'):
        label_str = ''
        het_SNVs = [(i.chrom, i.pos, i.ampl) for i in self.get_events('SNV')]
        red_SNVs = [] # wt -> het SNP -> LOH of ALT -> wt
        hom_SNVs = []
        for LOH_node in self.get_events('LOH'):
            for pos, allele in LOH_node.pos.items():
                if allele == 'REF':
                    hom_SNVs.append((LOH_node.chrom, pos))
                else:
                    if (LOH_node.chrom, pos, LOH_node.ampl) in het_SNVs:
                        n_ampl = sum(
                            [i.ampl == LOH_node.ampl for i in self.get_events() \
                                if i != LOH_node])
                        if n_ampl == 1:
                            red_SNVs.append((LOH_node.chrom, pos))

        for _, events in self.events.items():
            for event in sorted(events):
                if isinstance(event, SNV) \
                        and (event.chrom, event.pos) in hom_SNVs + red_SNVs:
                    continue
                elif isinstance(event, LOH) \
                        and (event.chrom, list(event.pos.keys())[0]) in red_SNVs:

                    continue
                label_str += f'{str(event)}{sep}'

        return f'{self.name} [label=<{label_str}>];'


class CellNode(Node):
    def __init__(self, name, style):
        super().__init__(name, style, False)
        self.cells = int(re.search('label="(\d+) cells', style).group(1))
        self.perc = int(re.search('\\\\n(\d+)\\\\%', style).group(1))
        self.cell_style = style.split('"')[-1].strip(' []')
        self.color = re.search('color=(\w+)', self.cell_style).group(1)
        self.cell_style = self.cell_style.replace(self.color, f'"{self.color}"')


    def get_percentage(self):
        return self.perc


    def get_cells(self):
        return self.cells


    def get_label(self):
        return f'{self.perc}\\%\\n{self.cells} cells'


    def add_cells(self, cells):
        self.cells += cells


    def add_percentage(self, perc):
        self.perc += perc


    def set_color(self, new_color):
        self.cell_style = self.cell_style.replace(self.color, new_color)
        self.color = new_color


    def __str__(self):
        return f'{self.name} [label="{self.get_label()}" {self.cell_style}];'


class Event:
    pass


class SNV(Event):
    def __init__(self, ampl, chrom, pos):
        self.ampl = ampl
        self.chrom = chrom
        self.pos = pos


    def get_pos_id(self):
        return (self.chrom, self.pos)


    def __lt__(self, other):
        # SNP Event
        try:
            return (CHROM[self.chrom], self.pos) \
                < (CHROM[other.chrom], other.pos)
        # LOH Event
        except TypeError:
            return (CHROM[self.chrom], self) \
                < (CHROM[other.chrom], list(other.pos.keys())[0])


    def __str__(self):
        # return f'{self.ampl}(chr{self.chrom}_{self.pos})'
        return f'chr{self.chrom}:{self.pos}'


class LOH(Event):
    def __init__(self, ampl, chrom, pos):
        self.ampl = ampl
        self.chrom = chrom
        self.pos = pos


    def __lt__(self, other):
        # SNP Event
        try:
            return (CHROM[self.chrom], list(self.pos.keys())[0]) \
                < (CHROM[other.chrom], other.pos)
        # LOH Event
        except TypeError:
            return (CHROM[self.chrom], list(self.pos.keys())[0]) \
                < (CHROM[other.chrom], list(other.pos.keys())[0])


    def __str__(self):
        # pos_str = ','.join(
        #     [f'chr{self.chrom}_{i[0]}-{i[1]}' for i in self.pos.items()])
        # return f'<B>LOH {self.ampl}({pos_str})</B>'
        pos_str = ''
        for i, j in self.pos.items():
            if j == 'ALT':
                pos_str += f',chr{self.chrom}:{i} (LOH ALT)'
            else:
                pos_str += f',chr{self.chrom}:{i} (HOM)'
        return pos_str.lstrip(',')


class Edge:
    def __init__(self, start, end, style):
        self.start = start
        self.end = end
        self.style = style.rstrip(';')


    def __str__(self):
        return f'{self.start.name} -> {self.end.name} {self.style};'


class CellEdge(Edge):
    def __init__(self, start, end, style):
        super().__init__(start, end, style)
        self.color = re.search('color=(\w+)', self.style).group(1)
        self.style = self.style.replace(self.color, f'"{self.color}"')

    def set_color(self, new_color):
        self.style = self.style.replace(self.color, new_color)
        self.color = new_color


class EventEdge(Edge):
    def __init__(self, start, end, style):
        super().__init__(start, end, style)


# ------------------------------------------------------------------------------


def main(args):
    for subdir, dirs, files in tqdm(os.walk(args.input)):
        for file in files:
            if not file.endswith('.gv'):
                continue
            file_abs = os.path.join(subdir, file)
            t = Tree(file_abs)
            t.remove_sparse_nodes(args.min_cells)
            if re.search('(?<=/)cl\d$|S\d$', subdir.rstrip('/')):
                cl_id = re.search('(?<=/)cl\d$|S\d$', subdir).group()
                t.set_cell_node_color(COLORS[cl_id])
            t.render(overwrite=args.overwrite)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        help='Input directory containing .gv files.')
    parser.add_argument('-mc', '--min_cells', type=float, default=1,
        help='Minimum percentage of cells for a node to be plotted. Default = 1.')
    parser.add_argument('-o', '--overwrite', action='store_true', default=False,
        help='If set, existing plots are overwritten.')

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args)