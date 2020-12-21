# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a Multilayer network (HS, MM & DM) and extracts only conserved genes.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')

    #
    #
    #
    print('Reading Network')
    rGfile_gpickle = 'results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Selecting conserved only
    #
    print('Selecting nodes with cross links to both other layers')
    print('Number of nodes: {:,d}'.format(G.number_of_nodes()))
    nodes_to_remove = []
    for node, nodedata in G.nodes(data=True):
        include = True
        set_connects_to_layers = set()
        for node, neighbor, edgedata in G.edges(node, data=True):
            if edgedata['type'] == 'cross':
                neighbor_layer = G.nodes[neighbor]['layer']
                set_connects_to_layers.add(neighbor_layer)
                if len(set_connects_to_layers) >= 2:
                    include = False
        if include:
            nodes_to_remove.append(node)
    G.remove_nodes_from(nodes_to_remove)
    #
    print('Number of nodes: {:,d}'.format(G.number_of_nodes()))

    # Removal of Nodes with only type=='cross' edges
    print('Removing nodes (with only cross edges)')
    print('Number of nodes: {:,d}'.format(G.number_of_nodes()))
    remove_isolates_nodes = []
    for node in G.nodes():
        if not any([True for i, j, d in G.edges(node, data=True) if d['type'] == 'intra']):
            remove_isolates_nodes.append(node)
    G.remove_nodes_from(remove_isolates_nodes)
    print('Number of nodes: {:,d}'.format(G.number_of_nodes()))

    ##
    # Export
    ##
    print('Exporting')
    wGfile_gpickle = 'results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
    ensurePathExists(wGfile_gpickle)
    nx.write_gpickle(G, wGfile_gpickle)
