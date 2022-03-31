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
    # Core Genes
    #
    rCfile = '../../02-core_genes/results/pipeline-core/{layer:s}_meiotic_genes.csv'.format(layer=layer)
    dfc = pd.read_csv(rCfile, index_col=0, usecols=['id_gene', 'gene'])
    #is_core = dfc.index.to_list()

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
    is_ortho_metric_edges = [(i, j) for i, j, d in OB.edges(data=True) if d.get('is_metric_ortho') == is_metric_ortho_string]
    is_ortho_metric_nodes = set(list(chain(*is_ortho_metric_edges)))

    #
    # Results
    #
    dfc['is_core'] = True
    dfc['is_conserved'] = dfc.index.isin(is_conserved)
    dfc['is_ortho_metric_nodes'] = dfc.index.isin(is_ortho_metric_nodes)
    #
    # Export
    dfc.to_csv('results/identity-{layer:s}-conserved-genes-core-ortho.csv'.format(layer=layer))
