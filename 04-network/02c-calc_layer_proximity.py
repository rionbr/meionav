# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and calculats a similarity across layers using gene homology.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
from itertools import combinations
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')

    print('Reading Network')
    rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    r = []
    for (layer_i, Gi), (layer_j, Gj) in combinations([('HS', HSG), ('MM', MMG), ('DM', DMG)], 2):
        print("Comparing: {layer_i:s} with {layer_j:s}".format(layer_i=layer_i, layer_j=layer_j))

        genes_i = [*Gi.nodes()]
        genes_j = [*Gj.nodes()]

        genes_ij = genes_i + genes_j

        # Only genes in this modules
        Gtmp = nx.subgraph(G, genes_ij).copy()

        # Remove intra edges
        remove_intra_edges = [(i, j) for i, j, d in Gtmp.edges(data=True) if d.get('type', None) == 'intra']
        Gtmp.remove_edges_from(remove_intra_edges)

        # Remove isolates
        remove_isolates_nodes = list(nx.isolates(Gtmp))
        Gtmp.remove_nodes_from(remove_isolates_nodes)

        a = set(genes_i)
        b = set(genes_j)
        a_union_b = a.union(b)
        a_inter_b = set(Gtmp.nodes())

        dist = len(a_inter_b) / len(a_union_b)
        r.append((layer_i, layer_j, dist))

    dfR = pd.DataFrame(r, columns=['layer-i', 'layer-j', 'proximity'])

    ##
    # Export
    ##
    print('Exporting')
    wCSVfile = 'results/proximity/proximity-{celltype:s}-{network:s}-{threshold:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str)
    ensurePathExists(wCSVfile)
    dfR.to_csv(wCSVfile)
