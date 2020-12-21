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

    data = {
        'HS': {
            'name': 'Human',
        },
        'MM': {
            'name': 'Mouse',
        },
        'DM': {
            'name': 'Insect',
        }
    }

    for celltype in celltypes:
        print("Loading celltype: {celltype:s}".format(celltype=celltype))
        #
        rGfile = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        G = nx.read_gpickle(rGfile)
        #
        for layer in layers:
            print('Separate layer {layer:s}'.format(layer=layer))
            Gl = get_network_layer(G, layer)
            #
            data[layer][celltype] = {
                'graph': Gl
            }

    # Compute Jaccard
    r = []
    for celltype_i, celltype_j in combinations(celltypes, 2):

        for layer in layers:

            G_i = data[layer][celltype_i]['graph']
            G_j = data[layer][celltype_j]['graph']

            genes_i = G_i.nodes()
            genes_j = G_j.nodes()

            # JAccard
            a = set(genes_i)
            b = set(genes_j)
            #
            a_union_b = a.union(b)
            a_inter_b = a.intersection(b)
            #
            prox = len(a_inter_b) / len(a_union_b)
            r.append((celltype_i, celltype_j, layer, prox))
    df = pd.DataFrame(r, columns=['celltype_i', 'celltype_j', 'layer', 'proximity'])

    # Export
    wCSVfile = 'results/network-jaccard-celltypes.csv'
    ensurePathExists(wCSVfile)
    df.to_csv(wCSVfile)
