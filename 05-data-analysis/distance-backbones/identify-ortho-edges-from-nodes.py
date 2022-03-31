# coding=utf-8
# Author: Rion B Correia
# Date: June 07, 2021
#
# Description: Identify conserved genes in the ortho-backbone
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
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    layer = 'DM'
    layers = ['HS', 'MM', 'DM']

    #
    # Conserved Network
    #
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)
    Gt = get_network_layer(G, layer)
    is_conserved = Gt.nodes()

    #
    # Ortho-backbone network
    #
    rOBfile_gpickle = '../../04-network/results/network-closure-ortho/{celltype:s}/net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network='thr', threshold=threshold_str, layer=layer)
    OB = nx.read_gpickle(rOBfile_gpickle)

    is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])
    #is_ortho_metric_edges = [(i, j) for i, j, d in OB.edges(data=True) if d.get('is_metric_ortho') == is_metric_ortho_string]
    #is_ortho_metric_nodes = set(list(chain(*is_ortho_metric_edges)))

    # FBgn0003943 = 
    nodes = [
        #'FBgn0003943',  # Ubi-p63E
        #'FBgn0014342',  # mia
        #'FBgn0002842',  # sa
        #'FBgn0002673',  # twe
        'FBgn0036516',  # CG7656
    ]
    
    for i in nodes:
        r = []
        for j, d in OB[i].items():
            gene_i = OB.nodes[i]['label']
            gene_j = OB.nodes[j]['label']
            is_conserved_i = i in is_conserved
            is_conserved_j = j in is_conserved
            weight_ij = d.get('weight')
            is_metric = d.get('is_metric')
            is_metric_ortho = d.get('is_metric_ortho')
            if is_metric:
                r.append((layer, i, j, gene_i, gene_j, is_conserved_i, is_conserved_j, weight_ij, is_metric, is_metric_ortho))
        # Results
        df = pd.DataFrame(r, columns=['layer', 'id-i', 'id-j', 'gene-i', 'gene-j', 'is_conserved-i', 'is_conserved-j', 'weight', 'is_metric', 'is_metric_ortho'])
        #Export
        df.to_csv('results/identify-orthoedges-from-genes-{layer:s}-{gene_i:s}.csv'.format(layer=layer, gene_i=gene_i))