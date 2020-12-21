# coding=utf-8
# Author: Rion B Correia
# Date: July 22, 2020
#
# Description: Reads threshold networks and calculates the jaccard measure on gene sets
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
from itertools import combinations


if __name__ == '__main__':

    celltypes = ['spermatocyte', 'enterocyte', 'muscle', 'neuron']
    #
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    layers = ['HS', 'MM', 'DM']

    data = {}

    for celltype in celltypes:
        print("Loading celltype: {celltype:s}".format(celltype=celltype))
        #
        rGfile = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        G = nx.read_gpickle(rGfile)
        #
        data[celltype] = {
            'graph': G
        }
        #
        for layer in layers:
            print('Separate layer {layer:s}'.format(layer=layer))
            Gl = get_network_layer(G, layer)
            #
            data[celltype][layer] = {
                'graph': Gl
            }

    # Compute Jaccard
    r = []
    for layer_i, layer_j in combinations(layers, 2):

        for celltype in celltypes:

            G = data[celltype]['graph']
            G_i = data[celltype][layer_i]['graph']
            G_j = data[celltype][layer_j]['graph']

            genes_i = list(G_i.nodes())
            genes_j = list(G_j.nodes())
            genes_ij = genes_i + genes_j
            #
            G_tmp = nx.subgraph(G, genes_ij).copy()

            # Remove intra edges
            remove_intra_edges = [(i, j) for i, j, d in G_tmp.edges(data=True) if d.get('type', None) == 'intra']
            G_tmp.remove_edges_from(remove_intra_edges)

            # Remove isolates
            remove_isolates_nodes = list(nx.isolates(G_tmp))
            G_tmp.remove_nodes_from(remove_isolates_nodes)

            # Jaccard
            a = set(genes_i)
            b = set(genes_j)
            #
            a_union_b = a.union(b)
            a_inter_b = set(G_tmp.nodes())
            #
            prox = len(a_inter_b) / len(a_union_b)
            r.append((layer_i, layer_j, celltype, prox))
    df = pd.DataFrame(r, columns=['layer_i', 'layer_j', 'celltype', 'proximity'])

    # Export
    wCSVfile = 'results/network-jaccard-species.csv'
    ensurePathExists(wCSVfile)
    df.to_csv(wCSVfile)
