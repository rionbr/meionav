# coding=utf-8
# Author: Rion B Correia
# Date: March 17, 2022
#
# Description: Computes the number of edges per attribute in conserved vs non-conserved nodes
#
#
import networkx as nx
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import get_network_layer, ensurePathExists


if __name__ == '__main__':

    celltype = 'spermatocyte'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    layers = ['HS', 'MM', 'DM']

    # Thresholded Multi-layer Network
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='thr', threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    # Conserved Multi-layer Network
    rCGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
    CG = nx.read_gpickle(rCGfile_gpickle)

    # identify conserved nodes
    is_conserved = set(CG.nodes())
    #dict_is_conserved = {n: True if n in is_conserved else False for n in G.nodes()}
    #nx.set_node_attributes(G, name='is_conserved', values=dict_is_conserved)

    #
    # Only for Human
    #
    r = []
    for layer in layers:
    

        Gt = get_network_layer(G, layer)

        df = nx.to_pandas_edgelist(Gt)
        df['i_conserved'] = df['source'].map(dict_is_conserved)
        df['j_conserved'] = df['target'].map(dict_is_conserved)

        cons_vs_cons = df.loc[(df['i_conserved'] == True) & (df['j_conserved'] == True), 'weight'].shape[0]
        cons_vs_nonc = df.loc[((df['i_conserved'] == True) & (df['j_conserved'] == False) | (df['i_conse rved'] == False) & (df['j_conserved'] == True)), 'weight'].shape[0]
        nonc_vs_nonc = df.loc[(df['i_conserved'] == False) & (df['j_conserved'] == False), 'weight'].shape[0]

        total = df.shape[0]
        print("Cons vs Cons: {n:d} ({p:.2%})".format(n=cons_vs_cons, p=cons_vs_cons / total))
        print("Cons vs Non-Cons: {n:d} ({p:.2%})".format(n=cons_vs_nonc, p=cons_vs_nonc / total))
        print("Non-Cons vs Non-Cons: {n:d} ({p:.2%})".format(n=nonc_vs_nonc, p=nonc_vs_nonc / total))

        r.append((layer, 'all-edges', cons_vs_cons, cons_vs_nonc, nonc_vs_nonc, total))

        for attribute in ['experiments', 'coexpression', 'database', 'textmining', 'cooccurence', 'fusion']:

            dft = df.loc[df[attribute].notnull(), :]
            #
            cons_vs_cons = dft.loc[(dft['i_conserved'] == True) & (dft['j_conserved'] == True), 'weight'].shape[0]
            cons_vs_nonc = dft.loc[((dft['i_conserved'] == True) & (dft['j_conserved'] == False) | (dft['i_conserved'] == False) & (dft['j_conserved'] == True)), 'weight'].shape[0]
            nonc_vs_nonc = dft.loc[(dft['i_conserved'] == False) & (dft['j_conserved'] == False), 'weight'].shape[0]

            total = dft.shape[0]

            r.append((layer, attribute, cons_vs_cons, cons_vs_nonc, nonc_vs_nonc, total))

    dfR = pd.DataFrame(r, columns=['layer', 'edge-type', 'cons-vs-cons', 'cons-vs-nonc', 'nonc-vs-nonc', 'total'])

    # Export
    dfR.to_csv('results/stats-edge-attribute-conserved.csv')
