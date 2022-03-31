# coding=utf-8
# Author: Rion B Correia
# Date: March 16, 2022
#
# Description: Calculates conserved gene homophily in networks and backbones
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
    layers = ['HS', 'MM', 'DM']
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    r = []
    for layer in layers:
        print("Layer:", layer)
        #
        # Threshold Network
        #
        rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='thr', threshold=threshold_str)
        G = nx.read_gpickle(rGfile_gpickle)
        Gt = get_network_layer(G, layer)

        n_nodes = Gt.number_of_nodes()
        n_edges = Gt.number_of_edges()
        # print("Thr:", n_nodes_thr, n_edges_thr)

        #
        # Conserved Network
        #
        rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
        CG = nx.read_gpickle(rGfile_gpickle)
        CGt = get_network_layer(CG, layer)
        is_conserved = set(CGt.nodes())

        n_nodes_cons = CGt.number_of_nodes()
        n_edges_cons = CGt.number_of_edges()

        dict_is_conserved = {n: True if n in is_conserved else False for n in Gt.nodes()}
        nx.set_node_attributes(Gt, name='is_conserved', values=dict_is_conserved)
        #
        att_assort_coef_thr = nx.attribute_assortativity_coefficient(Gt, attribute='is_conserved')

        #
        # Backbone
        #
        rBfile_gpickle = '../../04-network/results/network-closure/{celltype:s}/net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network='thr', threshold=threshold_str, layer=layer)
        B = nx.read_gpickle(rBfile_gpickle)
        not_metric = {(i, j) for i, j, d in B.edges(data=True) if d.get('is_metric') != True}
        B.remove_edges_from(not_metric)

        n_nodes_mb = B.number_of_nodes()
        n_edges_mb = B.number_of_edges()

        dict_is_conserved = {n: True if n in is_conserved else False for n in B.nodes()}
        nx.set_node_attributes(B, name='is_conserved', values=dict_is_conserved)
        #
        att_assort_coef_mb = nx.attribute_assortativity_coefficient(B, attribute='is_conserved')
        #
        # Output
        r.append((celltype, layer, n_nodes, n_edges, n_nodes_cons, n_edges_cons, n_nodes_mb, n_edges_mb, att_assort_coef_thr, att_assort_coef_mb))

    #
    columns = [
        'celltype',
        'layer',
        'n-nodes-thr',
        'n-edges-thr',
        'n-nodes-cons',
        'n-edges-cons',
        'n-nodes-metric',
        'n-edges-metric',
        'att-assort-coef',
        'att-assort-coef-backbone'
    ]
    df = pd.DataFrame(r, columns=columns)

    #Export
    df.to_csv('results/assortativity-conserved.csv')