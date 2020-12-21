# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and computed its backbone.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
# distanceclosure
from distanceclosure.dijkstra import Dijkstra
# from distanceclosure._dijkstra import _py_single_source_shortest_distances
from distanceclosure.cython._dijkstra import _cy_single_source_shortest_distances
from distanceclosure.utils import _prox2dist as prox2dist
import argparse


def compute_metric_backbone(G):
    print('> Computing Dijkstra APSP')
    dict_edges = {(i, j): prox2dist(d['weight']) for i, j, d in G.edges(data=True)}
    dij = Dijkstra.from_edgelist(dict_edges, directed=False, verbose=10)

    print('> Metric')
    poolresults = list(range(len(dij.N)))
    for node in dij.N:
        print('> Dijkstra node {node:d} of {n_nodes:d}'.format(node=node + 1, n_nodes=len(dij.N)))
        poolresults[node] = _cy_single_source_shortest_distances(node, dij.N, dij.E, dij.neighbors, ('min', 'sum'), verbose=2)
    shortest_distances, local_paths = map(list, zip(*poolresults))
    dij.shortest_distances = dict(zip(dij.N, shortest_distances))
    MSD = dij.get_shortest_distances(format='dict', translate=True)
    #
    dict_edges_backbone = {}
    dict_edges_s_values = {}
    for (i, j), d in dict_edges.items():

        # Backbone = distance == distance-closure
        if d == MSD[i][j]:
            dict_edges_backbone[(i, j)] = True
        else:
            dict_edges_backbone[(i, j)] = False

        # S-Value = distance / distance-closure
        dict_edges_s_values[(i, j)] = d / MSD[i][j]

    return dict_edges_backbone, dict_edges_s_values


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = 'thr'
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')

    #
    # Load Network
    #
    print('Reading {network:s} Network'.format(network=network))
    rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Selecting layers
    #
    print('Extracting SubGraphs')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    for layer, Gtmp in zip(['HS', 'MM', 'DM'], [HSG, MMG, DMG]):
        print('Computing Backbone ({layer:s})'.format(layer=layer))

        print('Create empty graph ({layer:s})'.format(layer=layer))
        B = nx.Graph()
        B.add_nodes_from(Gtmp.nodes())
        B.add_edges_from(Gtmp.edges())
        #
        # Compute Backbones
        #
        print('Dijkstra ({layer:s})'.format(layer=layer))
        dict_edges_backbone, dict_edges_s_values = compute_metric_backbone(Gtmp)

        # To DataFrame
        dfB = pd.DataFrame({'backbone': dict_edges_backbone, 's_values': dict_edges_s_values})

        ##
        # Export
        ##
        print('Exporting ({layer:s})'.format(layer=layer))
        wBfile = 'results/backbone/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-backbone.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        ensurePathExists(wBfile)
        dfB.to_csv(wBfile)
