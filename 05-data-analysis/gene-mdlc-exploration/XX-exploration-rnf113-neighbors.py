# coding=utf-8
# Author: Rion B Correia
# Date: Nov 24, 2021
#
# Description: 
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, ensurePathExists


if __name__ == '__main__':

    celltype = 'spermatocyte'
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    layer = 'DM'
    layers = ['HS', 'MM', 'DM']


    #
    # Network
    #
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    # Only orthorelations
    H = G.copy()
    edges2remove = [(i,j) for i,j,d in H.edges(data=True) if d.get('type', None) != 'cross']
    H.remove_edges_from(edges2remove)

    #
    # Add Backbone
    #
    path_backbone = '../../04-network/results/network-closure/{celltype:}/'.format(celltype=celltype)
    for layer in layers:
        rBfile_gpickle = path_backbone + 'net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        B = nx.read_gpickle(rBfile_gpickle)
        #
        is_metric = nx.get_edge_attributes(B, 'is_metric')
        #is_ultrametric = nx.get_edge_attributes(B, 'is_ultrametric')
        nx.set_edge_attributes(G, values=is_metric, name='is_metric')
        #nx.set_edge_attributes(G, values=is_ultrametric, name='is_ultrametric')

    #
    # Add Ortho-backbone
    #
    path_ortho_backbone = '../../04-network/results/network-closure-ortho/{celltype:s}/'.format(celltype=celltype)
    for layer in layers:
        rOfile_gpickle = path_ortho_backbone + 'net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        O = nx.read_gpickle(rOfile_gpickle)
        #
        is_metric_ortho = nx.get_edge_attributes(O, 'is_metric_ortho')
        nx.set_edge_attributes(G, values=is_metric_ortho, name='is_metric_ortho')

    # Remove non-metric edges
    edges2remove = [(i,j) for i, j, d in Gt.edges(data=True) if d.get('is_metric', None) != True]
    Gt.remove_edges_from(edges2remove)


    #layer = 'HS'
    #Gt = get_network_layer(G, layer)

    #search

    # HS
    # 'ENSG00000139797' = 'RNF113B'

    # MM
    #'ENSMUSG00000098134' = 'Rnf113a2'

    # DM
    # mdlc = 'FBgn0038772'

    for n, d in H.nodes(data=True):
        if d.get('label') == 'RNF-113':
            print(n, d)
            break

    i = 'ENSG00000139797'
    for j in nx.neighbors(Gt, i):
        di = Gt.nodes[i]
        dj = Gt.nodes[j]
        dij = Gt[i][j]
        print(i, di['label'], dj['label'])
        #print(dij)

    # Only neighbor genes to RNF113 (in the 3 species) 
    all_a = set()
    for i in ['ENSG00000139797', 'ENSMUSG00000098134', 'FBgn0038772']:
        a = set(nx.neighbors(G, i))
        all_a.update(a)
    X = nx.subgraph(G, nbunch=all_a).copy()
    
    # that also connect to another species
    nodes2remove = []
    for i in X.nodes():
        hascross = False
        for j in nx.neighbors(X, i):
            if X[i][j]['type'] == 'cross':
                hascross = True
        if not hascross:
            nodes2remove.append(i)

    X.remove_nodes_from(nodes2remove)

    # cross edges with low weight
    for i, j, d in X.edges(data=True):
        if d['type'] == 'cross':
            X[i][j]['weight'] = 0.20


    nx.write_graphml(X, 'test.graphml')