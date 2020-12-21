# coding=utf-8
# Author: Rion B Correia
# Date: July 22, 2020
#
# Description: Reads a threshold and a conserved multiLayer network to identify network patterns.
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
import random
from utils import ensurePathExists, get_network_layer
from itertools import combinations
from collections import defaultdict
import random


def random_partition_graph(blocks={'conserved': [1, 2], 'non-conserved': [3, 4]}, m=0, p=0.5):
    """ Generate a Random Partition Graph
    blocks = dict of blocks and their sizes
    m = number of edges
    p = probability of a node to connect to another node inside the same block
    """
    G = nx.Graph()
    # Add nodes
    node = 0
    block = 0
    nodes = defaultdict(list)
    for label, sizes in blocks.items():
        #
        for size in sizes:
            tmp_nodes = list(range(node, node + size))
            G.add_nodes_from(tmp_nodes, **{'block': block, 'type': label})
            nodes[block].extend(tmp_nodes)
            node += size
            block += 1
    # Add edges
    edgequeue = list(range(m))
    while len(edgequeue):
        # Intra
        if random.random() < p:
                bi = random.choice(list(nodes.keys()))
                ni, nj = random.sample(nodes[bi], k=2)
        # Cross
        else:
            bi, bj = random.sample(nodes.keys(), k=2)
            ni, nj = random.choice(nodes[bi]), random.choice(nodes[bj])
        #
        if not G.has_edge(ni, nj):
            G.add_edge(ni, nj)
            edgequeue.pop()
    #
    return G

def random_graph(n, m, blocks={'conserved': 1, 'non-conserved': 9}):
    """ Generate a Random Graph
    n = number of nodes
    m = number of edges
    block = number of nodes in the block
    """
    G = nx.gnm_random_graph(n=n, m=m)
    #
    assert sum(blocks.values()) == n , 'Number of nodes must be qual to the sum of the block sizes.'
    #
    nodes = set(list(G.nodes()))
    for label, size in blocks.items():

        block_nodes = random.sample(nodes, k=size)
        # Remove
        nodes = nodes.difference(block_nodes)
        #
        values = {node: label for node in block_nodes}
        nx.set_node_attributes(G, values=values, name='type')
    return G


def node_attribute_assortativity(G, attribute):
    r = {}
    for node, data in G.nodes(data=True):
        node_attribute = data[attribute]
        neighbors = list(G.neighbors(node))
        equals = len([neighbor for neighbor in neighbors if G.nodes[neighbor][attribute] == node_attribute])
        totals = len(neighbors)
        assortativity = equals / totals
        r[node] = assortativity
    return r



if __name__ == '__main__':

    celltype = 'spermatocyte'

    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    layers = ['HS', 'MM', 'DM']

    rGtfile = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='thr', threshold=threshold_str)
    Gt = nx.read_gpickle(rGtfile)

    rGcfile = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
    Gc = nx.read_gpickle(rGcfile)

    # Remove Cross edges
    Gt.remove_edges_from([(i, j) for i, j, d in Gt.edges(data=True) if d['type'] == 'cross'])

    # Identify conserved nodes
    values = {k: 'conserved' if k in Gc.nodes() else 'non-conserved' for k in Gt.nodes()}
    nx.set_node_attributes(Gt, values=values, name='type')

    data = {
        'HS': {
            'name': 'Human',
            'facecolor': '#2ca02c',
            'edgecolor': '#98df8a',
        },
        'MM': {
            'name': 'Mouse',
            'facecolor': '#7f7f7f',
            'edgecolor': '#c7c7c7',
        },
        'DM': {
            'name': 'Insect',
            'facecolor': '#ff7f0e',
            'edgecolor': '#ffbb78',
        }
    }

    r = []
    nr = {}
    for layer in layers:
        print('Layer: {layer:s}'.format(layer=layer))
        Gtl = get_network_layer(Gt, layer)
        Gcl = get_network_layer(Gc, layer)
        density = nx.density(Gtl)
        assort = nx.attribute_assortativity_coefficient(Gtl, attribute='type')
        #
        r.append((layer, None, density, None, assort))
        #
        node_assort = node_attribute_assortativity(Gtl, attribute='type')
        nr[('DM', None)] = list(node_assort.values())
        #
        data[layer]['graph-thr'] = Gtl
        data[layer]['graph-conserved'] = Gcl
        data[layer]['assortativity'] = assort

    #
    for layer in layers:
        Gtl = data[layer]['graph-thr']
        Gcl = data[layer]['graph-conserved']
        #
        number_of_nodes = Gtl.number_of_nodes()
        number_of_edges = Gtl.number_of_edges()
        number_of_nodes_conserved = Gcl.number_of_nodes()
        
        assort = nx.attribute_assortativity_coefficient(Gtl, attribute='type')
        
        md = nx.attribute_mixing_dict(Gtl, attribute='type', normalized=False)
        r.append((layer, assort, md['conserved']['conserved'], md['conserved']['non-conserved'], md['non-conserved']['conserved'], md['non-conserved']['non-conserved']))
        #$\
    df = pd.DataFrame(r, columns=['layer', 'assortativity', 'con-con', 'con-non', 'non-con', 'non-non'])
    print(df)

        # RANDOM NULL MODEL
        """
        for i in range(5):
            Gtmp = random_graph(number_of_nodes, number_of_edges, {'conserved': number_of_nodes_conserved, 'non-conserved': (number_of_nodes - number_of_nodes_conserved)})
            node_assort = node_attribute_assortativity(Gtmp, attribute='type')
            nr[('DM', i)] = list(node_assort.values())
        dfnr = pd.DataFrame(nr)

        wCSVFile = 'results/csv-null-model-node-assortativity.csv'
        ensurePathExists(wCSVFile)
        dfnr.to_csv(wCSVFile)
        """

        """
        # RANDOM PARTITION NULL MODEL
        number_of_blocks = 50
        search_space = np.linspace(0, 1, 25, endpoint=True)
        for p in search_space:
            print('p: {p:.6f}'.format(p=p))

            for i in range(5):
                print('iteration: {i:d}'.format(i=i))
                #n_m_major = math.floor(number_of_nodes * 0.50 / (0.20 * b))
                #_m_minor = math.floor(number_of_nodes * 0.50 / (0.80 * b))
                #sizes = [n_m_major] * math.floor(.20 * b) + [n_m_minor] * math.floor(.80 * b)
                #
                number_of_nodes_per_block = math.floor(number_of_nodes / number_of_blocks)
                block_sizes = [number_of_nodes_per_block] * number_of_blocks
                # add the difference in the last block
                diff_number_of_nodes = number_of_nodes - sum(block_sizes)
                block_sizes[-1] += diff_number_of_nodes

                # Select random partitions to be set as conserved
                tmp_blocks_sizes = list(block_sizes)
                tmp_blocks_sizes_conserved = list()
                #random.shuffle(m_ids)
                while sum(tmp_blocks_sizes_conserved) < number_of_nodes_conserved:
                    tmp_blocks_sizes_conserved.append(tmp_blocks_sizes.pop(0))

                tmp_number_of_blocks_conserved = len(tmp_blocks_sizes_conserved)
                tmp_number_of_blocks_non_conserved = len(tmp_blocks_sizes)
                #
                tmp_blocks = {'conserved': tmp_blocks_sizes_conserved, 'non-conserved': tmp_blocks_sizes}
                #Gtmp = nx.random_partition_graph(block_sizes, p_in=p_in_adj, p_out=p_out_adj, directed=False)
                Gtmp = random_partition_graph(tmp_blocks, m=number_of_edges, p=p)
                tmp_density = nx.density(Gtmp)
                #
                assort = nx.attribute_assortativity_coefficient(Gtmp, attribute='type')
                #
                r.append((layer, i, tmp_density, p, assort))
        """
    dfr = pd.DataFrame(r, columns=['layer', 'i', 'density', 'p', 'assortativity'])
    print(dfr)

    # Export
    wCSVFile = 'results/csv-null-model-assortativity.csv'
    ensurePathExists(wCSVFile)
    df.to_csv(wCSVFile)
