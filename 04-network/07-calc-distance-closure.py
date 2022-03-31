# coding=utf-8
# Author: Rion B Correia
# Date: March 23, 2021
#
# Description: Reads a MultiLayer network (HS, MM & DM) and computed its modules using Louvain & Infomap.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
# Distance Closure
from distanceclosure import all_pairs_dijkstra_path_length
from distanceclosure.cython.dijkstra import cy_all_pairs_dijkstra_path_length
from distanceclosure.utils import prox2dist, from_networkx_to_dijkstra_format
#
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    #parser.add_argument("--layer", default='DM', type=str, choices=['DM', 'MM', 'HS'], help="Network layer to compute SVD. Defaults to 'DM'.")

    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    #layer = args.layer

    #
    # Load Network
    #
    print('Reading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
    rGfile_gpickle = 'results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    # SVD per Layer
    for layer in ['HS', 'MM', 'DM']:
        print('Isolate {layer:s} Layer'.format(layer=layer))
        Gt = get_network_layer(G, layer=layer)
        #
        # Prox to Distance
        P_dict = nx.get_edge_attributes(Gt, 'weight')
        D_dict = {key: prox2dist(value) for key, value in P_dict.items()}
        nx.set_edge_attributes(Gt, name='distance', values=D_dict)

        # Dijkstra Format
        nodes, edges, neighbors = from_networkx_to_dijkstra_format(Gt, weight='distance')
        dict_int_nodes = {i: u for i, u in enumerate(Gt.nodes())}
        #
        # Metric computation
        #
        c = 1
        #for i, metric_distances in all_pairs_dijkstra_path_length(Gt, weight='distance', disjunction=sum):
        for i, metric_distances in cy_all_pairs_dijkstra_path_length(nodes, edges, neighbors, ('min', 'sum')):
            i = dict_int_nodes[i]
            metric_distances = {dict_int_nodes[j]: d for j, d in metric_distances.items()}

            if c % 10 == 0:
                print('> {layer:s} Metric Dijkstra: {c:} of {total:}'.format(layer=layer, c=c, total=Gt.number_of_nodes()))
            for j, metric_distance in metric_distances.items():

                # New Edge?
                if not Gt.has_edge(i, j):
                    # Self-loops have proximity 1, non-existent have 0
                    """
                    proximity = 1.0 if i == j else 0.0
                    Gt.add_edge(i, j, distance=np.inf, proximity=proximity, metric_distance=float(cm), is_metric=False)
                    """
                else:
                    Gt[i][j]['metric_distance'] = metric_distance
                    Gt[i][j]['is_metric'] = True if ((metric_distance == Gt[i][j]['distance']) and (metric_distance != np.inf)) else False
            c += 1

        #
        # Ultrametric computation
        #
        c = 1
        #for i, ultrametric_distances in all_pairs_dijkstra_path_length(Gt, weight='distance', disjunction=max):
        for i, ultrametric_distances in cy_all_pairs_dijkstra_path_length(nodes, edges, neighbors, ('min', 'max')):
            i = dict_int_nodes[i]
            ultrametric_distances = {dict_int_nodes[j]: d for j, d in ultrametric_distances.items()}

            if c % 10 == 0:
                print('> {layer:s} Ultrametric Dijkstra: {c:} of {total:}'.format(layer=layer, c=c, total=Gt.number_of_nodes()))
            for j, ultrametric_distance in ultrametric_distances.items():

                # New Edge?
                if not Gt.has_edge(i, j):
                    # Self-loops have proximity 1, non-existent have 0
                    """
                    proximity = 1.0 if i == j else 0.0
                    Gt.add_edge(i, j, distance=np.inf, proximity=proximity, ultrametric_distance=float(cum), is_ultrametric=False)
                    """
                else:
                    Gt[i][j]['ultrametric_distance'] = ultrametric_distance
                    Gt[i][j]['is_ultrametric'] = True if ((ultrametric_distance == Gt[i][j]['distance']) and (ultrametric_distance != np.inf)) else False
            c += 1

        #
        # S-Value
        #
        print('--- Calculating S Values ---')

        S_dict = {
            (i, j): float(d['distance'] / d['metric_distance'])
            for i, j, d in Gt.edges(data=True)
            if ((d.get('distance') < np.inf) and (d.get('metric_distance') > 0))
        }
        nx.set_edge_attributes(Gt, name='s-value', values=S_dict)



        print('Saving results to .CSV')
        wGpickle = 'results/network-closure/{celltype:s}/net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        print('> gPickle')
        nx.write_gpickle(Gt, wGpickle)
        #
    print('Done.')
