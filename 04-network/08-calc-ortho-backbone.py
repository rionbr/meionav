# coding=utf-8
# Author: Rion B Correia
# Date: March 23, 2021
#
# Description: Reads a MultiLayer network (HS, MM & DM), its backbone, and computes the ortho-backbone.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from itertools import product
from utils import ensurePathExists, get_network_layer
#
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte'] #, 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")

    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    #
    layers = ['HS', 'MM', 'DM']
    #
    # Load Network
    #
    data = {
        'closure': {
            'HS': None,
            'MM': None,
            'DM': None
        },
        'backbone': {
            'HS': None,
            'MM': None,
            'DM': None
        }
    }
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    path_net = 'results/network/{celltype:}/'.format(celltype=celltype)
    rGfile = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile)

    is_not_cross_edges = [(i, j) for i, j, d in G.edges(data=True) if d.get('type') != 'cross']
    H = G.copy()
    H.remove_edges_from(is_not_cross_edges)

    #
    # Backbones
    #
    for layer in layers:
        path_backbone = 'results/network-closure/{celltype:}/'.format(celltype=celltype)
        rBfile = path_backbone + 'net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        C = nx.read_gpickle(rBfile)
        is_metric_edges = [(i, j) for i, j, d in C.edges(data=True) if d.get('is_metric') == True]
        B = C.edge_subgraph(is_metric_edges).copy()
        #
        data['closure'][layer] = C
        data['backbone'][layer] = B

    #
    # "Ortho metric backbone" (metric backbone edges present in another layer)
    #
    r = []
    for layer in layers:
        print("- Computing {layer:s} ortho backbone".format(layer=layer))
        C = data['closure'][layer]
        B = data['backbone'][layer]

        number_of_metric_edges = B.number_of_edges()

        dict_is_metric_ortho = {}
        # Loop all backbone edges
        for (i, j) in B.edges():

            at_least_one_metric_ortho = False
            is_metric_ortho_string = 'is_metric_ortho'
            for other_layer in layers:

                if other_layer == layer:
                    continue

                ortho_is = [u for u in H.neighbors(i) if H.nodes[u]['layer'] == other_layer]
                ortho_js = [u for u in H.neighbors(j) if H.nodes[u]['layer'] == other_layer]

                combinations = list(product(ortho_is, ortho_js))

                # at least one of the nodes does not have orthologs
                if len(combinations) == 0:
                    continue

                BO = data['backbone'][other_layer]

                is_metric_ortho = any([BO.has_edge(ortho_i, ortho_j) for ortho_i, ortho_j in combinations])

                if is_metric_ortho:
                    at_least_one_metric_ortho = True
                    is_metric_ortho_string += '-{other_layer:s}'.format(other_layer=other_layer)

            if at_least_one_metric_ortho:
                dict_is_metric_ortho[(i, j)] = is_metric_ortho_string
            else:
                dict_is_metric_ortho[(i, j)] = False

        #
        nx.set_edge_attributes(C, values=dict_is_metric_ortho, name='is_metric_ortho')
        #nx.set_edge_attributes(B, values=dict_is_metric_ortho, name='is_metric_ortho')

        wCpickle = 'results/network-closure-ortho/{celltype:s}/net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        print('> gPickle')
        ensurePathExists(wCpickle)
        nx.write_gpickle(C, wCpickle)
        #
    print('Done.')
