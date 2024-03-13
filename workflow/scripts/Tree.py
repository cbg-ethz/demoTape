#!/usr/bin/env python3

import copy
import re


CHROM = {**{str(i): i for i in range(1, 23, 1)}, **{'X': 23, 'Y': 24}}


class Tree:
    STYLE = '[color=dimgray fontsize=24 fontcolor=black fontname=Helvetica penwidth=5]'
    GV = 'digraph G{{\nnode {style};\n{mut_edges}\n{mut_nodes}\n{cell_edges}\n{cell_nodes}\n}}'


    def __init__(self):
        self.nodes = {}
        self.edges = []
        self.file = ''


    def read(self, tree_file):
        self.file = tree_file
        with open(tree_file, 'r') as f:
            tree_str = f.read().strip()
        # Init nodes
        for line in tree_str.split('\n'):
            if line.startswith(('digraph', 'node', '}')) or '->' in line:
                continue
            self.add_node(line)

        # Init edges
        for line in tree_str.split('\n'):
            if line.startswith(('digraph', 'node', '}')) or not '->' in line:
                continue
            self.add_edge(line)

        self.get_root() # Sanity check for 1 root node


    def add_node(self, line):
        break_idx = line.index('[')
        name = int(line[:break_idx].strip())
        style = line[break_idx:].rstrip(';')

        if 'style' in style:
            self.nodes[name] = CellNode(name, style)
        else:
            self.nodes[name] = EventNode(name, style)


    def copy_node(self, node):
        if 'style' in node.style:
            self.nodes[node.name] = CellNode(node.name, node.style)
        else:
            self.nodes[node.name] = EventNode(node.name, node.style)


    def add_edge(self, line):
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


    def copy_edge(self, edge):
        start = self.nodes[edge.start.name]
        end = self.nodes[edge.end.name]
        if isinstance(end, EventNode):
            end.not_root()

        if 'dir=none' in edge.style:
            self.edges.append(CellEdge(start, end, edge.style))
        else:
            self.edges.append(EventEdge(start, end, edge.style))


    def combine(self, trees):
        new_root_name = 0
        new_root = EventNode(new_root_name, '[label=<<br/>>]')
        self.nodes[new_root_name] = new_root

        for tree_id, tree in trees:
            for node in tree.get_nodes():
                self.copy_node(node)
                self.nodes[node.name].set_id(tree_id)
            for edge in tree.get_edges():
                self.copy_edge(edge)

            old_root = self.nodes[tree.get_root().name]
            old_root.not_root() # anymore
            self.edges.append(
                EventEdge(new_root, old_root, '[color=dimgray penwidth=4 weight=2]'))

        #     try:
        #         root_snvs = root_snvs.intersection(old_root.get_events('SNV''))
        #     except NameError:
        #         root_snvs = set(old_root.get_events('SNV'))

        # new_root.add_SNVs(root_snvs)
        # for _, tree in trees:
        #     old_root = self.nodes[tree.get_root().name]
        #     old_root.remove_SNVs(root_snvs)


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
                    if self.is_tip(edge.start):
                        self.nodes.pop(edge.start.name)
                        self.nodes.pop(edge.end.name)
                        del_edges.append(edge)
                        del_edges.append(self.get_upstream_edge(edge.start))

            for del_edge in del_edges:
                self.edges.remove(del_edge)

            self._adjust_indices()

            if not del_edges:
                break


    def no_nodes(self, n_type='all'):
        if n_type == 'all':
            return len(self.nodes)
        elif n_type == 'cell':
            return len([i for i in self.nodes.values() if isinstance(i, CellNode)])
        else:
            return len([i for i in self.nodes.values() if isinstance(i, EventNode)])


    def increase_node_id(self, mut_add, cell_add):
        new_nodes = {}
        for node_id, node in self.nodes.copy().items():
            old_id = node_id
            if isinstance(node, EventNode):
                new_id = int(node_id) + mut_add
            else:
                new_id = int(node_id) + cell_add
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


    def is_tip(self, node):
        is_tip = True
        for edge in self.get_edges('event'):
            if edge.start == node:
                is_tip = False
        return is_tip


    def to_gv(self, out_file):
        with open(out_file, 'w') as f:
            f.write(str(self))


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
        for _, events in self.events.items():
            for event in events:
                label_str += f'{str(event)}{sep}'
        return f'{self.name} [label=<{label_str}>];'


class CellNode(Node):
    def __init__(self, name, style):
        super().__init__(name, style, False)
        self.cells = int(re.search('label="(\d+) cells', style).group(1))
        self.perc = int(re.search('\\\\n(\d+)\\\\%', style).group(1))
        self.cell_style = style.split('"')[-1].strip(' []')


    def get_percentage(self):
        return self.perc


    def get_label(self):
        return f'{self.cells} cells\\n{self.perc}\\%'


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


    def __str__(self):
        return f'{self.ampl}(chr{self.chrom}_{self.pos})'


class LOH(Event):
    def __init__(self, ampl, chrom, pos):
        self.ampl = ampl
        self.chrom = chrom
        self.pos = pos


    def __str__(self):
        pos_str = ','.join(
            [f'chr{self.chrom}_{i[0]}-{i[1]}' for i in self.pos.items()])
        return f'<B>LOH {self.ampl}({pos_str})</B>'


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


class EventEdge(Edge):
    def __init__(self, start, end, style):
        super().__init__(start, end, style)


if __name__ == '__main__':
    print('Here be dragons....')