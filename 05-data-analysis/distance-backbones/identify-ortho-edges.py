# coding=utf-8
# Author: Rion B Correia
# Date: March 30, 2021
#
# Description: Identify ortho edges in the ortho-backbone
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from itertools import product, chain
from collections import Counter
from utils import get_network_layer, ensurePathExists

if __name__ == '__main__':

    celltype = 'spermatocyte'
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    layer = 'DM'
    layers = ['HS', 'MM', 'DM']
    #

    # Identify via file
    dfs = pd.read_csv('results/ortho-backbone-subgraph/net-ortho-backbone-HS-subgraph-GENES-etal.csv', index_col=0)
    edges = dfs[['id_gene_i', 'id_gene_j']].to_records(index=False).tolist()


    # Identify via selected edges
    """
    edge_labels = [
        ('CDYL', 'MIER1'),
        ('DAZAP1', 'LRRC34'),
        ('DDX25', 'GLE1'),
        ('SPATA5', 'AFG1L'),
        ('RSBN1', 'SPATA20'),
        ('RSBN1', 'SPAG7'),
        ('SPATA20', 'SPAG7'),
        #
        ('BOLL', 'CPEB1'),
        ('TEX2', 'DNAJC14'),
        ('TEX2', 'ARV1'),
        ('TEX2', 'NKIRAS1'),
        ('TEX2', 'AC023055.1'),
    ]
    #edges = []
    for label_i, label_j in edge_labels:
        id_i = [n for n, d in data['graphs']['HS'].nodes(data=True) if d.get('label', None) == label_i][0]
        id_j = [n for n, d in data['graphs']['HS'].nodes(data=True) if d.get('label', None) == label_j][0]
        edges.append((id_i, id_j))
    """
    
    #
    data = {
        'graphs': {
            'HS': None,
            'MM': None,
            'DM': None
        },
        'backbones': {
            'HS': None,
            'MM': None,
            'DM': None
        }
    }
    #
    # Network
    #
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    # Intra layers
    is_not_cross_edges = [(i, j) for i, j, d in G.edges(data=True) if d.get('type') != 'cross']
    H = G.copy()
    H.remove_edges_from(is_not_cross_edges)

    #
    path_backbone = '../../04-network/results/network-closure/{celltype:}/'.format(celltype=celltype)
    for layer in ['HS', 'MM', 'DM']:
        # Layer
        Gt = get_network_layer(G, layer)
        data['graphs'][layer] = Gt
        #
        # Backbone
        #
        rBfile_gpickle = path_backbone + 'net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        B = nx.read_gpickle(rBfile_gpickle)
        data['backbones'][layer] = B

    #
    # Identify ortho edges
    #
    r = []
    #layer = 'HS'
    print("- Computing {layer:s} ortho backbone".format(layer=layer))
    Gt = data['graphs'][layer]
    B = data['backbones'][layer]

    # Loop all backbone edges
    for (i, j) in edges:
        label_i = Gt.nodes[i]['label']
        label_j = Gt.nodes[j]['label']
        weight_ij = Gt[i][j]['weight']
        #
        for other_layer in layers:

            if other_layer == layer:
                continue

            ortho_is = [u for u in H.neighbors(i) if H.nodes[u]['layer'] == other_layer]
            ortho_js = [u for u in H.neighbors(j) if H.nodes[u]['layer'] == other_layer]

            combinations = list(product(ortho_is, ortho_js))

            # at least one of the nodes does not have orthologs
            if len(combinations) == 0:
                continue

            BO = data['backbones'][other_layer]

            for ortho_i, ortho_j in combinations:
                if BO.has_edge(ortho_i, ortho_j):
                    ortho_i_label = BO.nodes[ortho_i]['label']
                    ortho_j_label = BO.nodes[ortho_j]['label']
                    ortho_weight_ij = BO[ortho_i][ortho_j]['weight']
                    r.append((layer, i, j, label_i, label_j, weight_ij, other_layer, ortho_i, ortho_j, ortho_i_label, ortho_j_label, ortho_weight_ij))

    dfR = pd.DataFrame(r, columns=['layer', 'id-i', 'id-j', 'gene-i', 'gene-j', 'weight-ij', 'other-layer', 'id-ortho-i', 'id-ortho-j', 'gene-ortho-i', 'gene-ortho-j', 'ortho-weight-ij'])
    # Export
    dfR.to_csv('results/selected-ortho-edges-GENES.csv')
