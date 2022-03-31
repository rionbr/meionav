# coding=utf-8
# Author: Rion B Correia
# Date: Jan 01, 2022
#
# Description: What happened to the ortho-backbone afther the downregulation of mdlc
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, ensurePathExists
import random


if __name__ == '__main__':

    celltype = 'spermatocyte'
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    layer = 'HS'
    layers = ['HS', 'MM', 'DM']

    #
    # Network
    #
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)
    G = get_network_layer(G, layer)

    #
    # Add Backbone
    #
    path_backbone = '../../04-network/results/network-closure/{celltype:}/'.format(celltype=celltype)
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
    rOfile_gpickle = path_ortho_backbone + 'net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    O = nx.read_gpickle(rOfile_gpickle)
    #
    is_metric_ortho = nx.get_edge_attributes(O, 'is_metric_ortho')
    nx.set_edge_attributes(G, values=is_metric_ortho, name='is_metric_ortho')


    # DM
    # mdlc / prp19 DGE Results
    #
    if layer == 'DM':
        rMDLCFile = '../../01-diff-gene-exp/results/mdlc/{layer:s}-DGE-mdlc_vs_control.csv'.format(layer=layer)
        dfM = pd.read_csv(rMDLCFile, index_col=0, usecols=['id', 'gene', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'])
        # Filter only DGE significant
        dfM = dfM.loc[(dfM['logFC'].abs() > 1) & (dfM['FDR'] <= 0.05) & (dfM['logCPM'] >= 1), :].copy()
        dfM['up/down'] = dfM['logFC'].map(lambda x: 'up' if x > 0 else 'down')
        #
        is_mdlc_dge = dfM['up/down'].to_dict()
        #
        nx.set_node_attributes(G, name='is_mdlc_dge', values=is_mdlc_dge)

        rPRP19File = '../../01-diff-gene-exp/results/prp19/{layer:s}-DGE-prp19_vs_control.csv'.format(layer=layer)
        dfP = pd.read_csv(rPRP19File, index_col=0, usecols=['id', 'gene', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'])
        # Filter only DGE significant
        dfP = dfP.loc[(dfP['logFC'].abs() > 1) & (dfP['FDR'] <= 0.05) & (dfP['logCPM'] >= 1), :].copy()
        dfP['up/down'] = dfP['logFC'].map(lambda x: 'up' if x > 0 else 'down')
        #
        is_prp19_dge = dfP['up/down'].to_dict()
        #
        nx.set_node_attributes(G, name='is_prp19_dge', values=is_prp19_dge)

    #
    # HS
    # RNF113B DGE Results
    #
    if layer == 'HS':
        rRNFFile = '../../01-diff-gene-exp/results/rnf113/{layer:s}-DGE-rnf113_vs_control-SELECTED.csv'.format(layer=layer)
        dfM = pd.read_csv(rRNFFile, index_col=0, usecols=['id', 'gene', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'])
        # No need to filter as this file only contains the SELECTED genes

        dfM['up/down'] = dfM['logFC'].map(lambda x: 'up' if x > 0 else 'down')
        #
        is_rnf113_dge = dfM['up/down'].to_dict()
        nx.set_node_attributes(G, name='is_rnf113_dge', values=is_rnf113_dge)

    is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])

    r = []
    for i, j, d in G.edges(data=True):

        i_gene = G.nodes[i].get('label', None)
        j_gene = G.nodes[j].get('label', None)

        e_is_metric = d.get('is_metric', None)
        e_is_ortho = d.get('is_metric_ortho', None)

        if layer == 'DM':
            i_is_mdlc = G.nodes[i].get('is_mdlc_dge', None)
            j_is_mdlc = G.nodes[j].get('is_mdlc_dge', None)
            i_is_prp19 = G.nodes[i].get('is_prp19_dge', None)
            j_is_prp19 = G.nodes[j].get('is_prp19_dge', None)

            r.append((i, j, i_gene, j_gene, i_is_mdlc, j_is_mdlc, i_is_prp19, j_is_prp19, e_is_metric, e_is_ortho))
        elif layer == 'HS':
            i_is_rnf113 = G.nodes[i].get('is_rnf113_dge', None)
            j_is_rnf113 = G.nodes[j].get('is_rnf113_dge', None)

            r.append((i, j, i_gene, j_gene, i_is_rnf113, j_is_rnf113, e_is_metric, e_is_ortho))

    if layer == 'DM':
        df = pd.DataFrame(r, columns=['i', 'j', 'i-gene', 'j-gene', 'i-is_mdlc', 'j-is_mdlc', 'i-is_prp19', 'j-is_prp19', 'is_metric', 'is_ortho'])
    elif layer == 'HS':
        df = pd.DataFrame(r, columns=['i', 'j', 'i-gene', 'j-gene', 'i-is_rnf113', 'j-is_rnf113', 'is_metric', 'is_ortho'])

    #Export
    if layer == 'DM':
        df.to_csv("results/net-DM-backbone-after-mdlc-prp19.csv")
    elif layer == 'HS':
        df.to_csv("results/net-HS-backbone-after-rnf113.csv")

    # Edges
    n_edges = G.number_of_edges()
    # Edges metric
    n_is_metric = df.loc[(df['is_metric'] == True), :].shape[0]
    # Edges ortho-backbone
    n_is_ortho = df.loc[(df['is_ortho'] == is_metric_ortho_string), :].shape[0]
    print('n_edges:', n_edges)
    print('n_is_metric:', n_is_metric)
    print('n_is_ortho:', n_is_ortho)

    if layer == 'DM':
        n_is_mdlc_dge = len(is_mdlc_dge)
        n_is_prp19_dge = len(is_prp19_dge)
    elif layer == 'HS':
        n_is_rnf113_dge = len(is_rnf113_dge)

    # Count mdlc disrupted
    if layer == 'DM':
        n_d_mdlc = df.loc[((df['i-is_mdlc'].notnull()) | (df['j-is_mdlc'].notnull())), :].shape[0]
        n_d_mdlc_metric = df.loc[(df['is_metric'] == True) & ((df['i-is_mdlc'].notnull()) | (df['j-is_mdlc'].notnull())), :].shape[0]
        n_d_mdlc_ortho = df.loc[(df['is_ortho'] == is_metric_ortho_string) & ((df['i-is_mdlc'].notnull()) | (df['j-is_mdlc'].notnull())), :].shape[0]
        print('n_d_mdlc:', n_d_mdlc)
        print('n_d_mdlc_metric:', n_d_mdlc_metric)
        print('n_d_mdlc_ortho:', n_d_mdlc_ortho)

    # Count prp19 disrupted
        n_d_prp19 = df.loc[((df['i-is_prp19'].notnull()) | (df['j-is_prp19'].notnull())), :].shape[0]
        n_d_prp19_metric = df.loc[(df['is_metric'] == True) & ((df['i-is_prp19'].notnull()) | (df['j-is_prp19'].notnull())), :].shape[0]
        n_d_prp19_ortho = df.loc[(df['is_ortho'] == is_metric_ortho_string) & ((df['i-is_prp19'].notnull()) | (df['j-is_prp19'].notnull())), :].shape[0]
        print('n_d_prp19:', n_d_prp19)
        print('n_d_prp19_metric:', n_d_prp19_metric)
        print('n_d_prp19_ortho:', n_d_prp19_ortho)

    # Count rnf113 disrupted
    elif layer == 'HS':
        n_d_rnf113 = df.loc[((df['i-is_rnf113'].notnull()) | (df['j-is_rnf113'].notnull())), :].shape[0]
        n_d_rnf113_metric = df.loc[(df['is_metric'] == True) & ((df['i-is_rnf113'].notnull()) | (df['j-is_rnf113'].notnull())), :].shape[0]
        n_d_rnf113_ortho = df.loc[(df['is_ortho'] == is_metric_ortho_string) & ((df['i-is_rnf113'].notnull()) | (df['j-is_rnf113'].notnull())), :].shape[0]
        print('n_d_rnf113:', n_d_rnf113)
        print('n_d_rnf113_metric:', n_d_rnf113_metric)
        print('n_d_rnf113_ortho:', n_d_rnf113_ortho)

    # random samples
    gene_ids = G.nodes()
    #
    r = []
    for run in range(100):
        #
        if layer == 'DM':
            random_mdlc_gene_ids = set(random.sample(gene_ids, k=n_is_mdlc_dge))
            random_prp19_gene_ids = set(random.sample(gene_ids, k=n_is_prp19_dge))

            df['i-is_mdlc_random'] = df['i'].map(lambda x: x in random_mdlc_gene_ids)
            df['j-is_mdlc_random'] = df['j'].map(lambda x: x in random_mdlc_gene_ids)

            df['i-is_prp19_random'] = df['i'].map(lambda x: x in random_prp19_gene_ids)
            df['j-is_prp19_random'] = df['j'].map(lambda x: x in random_prp19_gene_ids)

            n_d_mdlc = df.loc[((df['i-is_mdlc_random'] == True) | (df['j-is_mdlc_random'] == True)), :].shape[0]
            n_d_mdlc_metric = df.loc[(df['is_metric'] == True) & ((df['i-is_mdlc_random'] == True) | (df['j-is_mdlc_random'] == True)), :].shape[0]
            n_d_mdlc_ortho = df.loc[(df['is_ortho'] == is_metric_ortho_string) & ((df['i-is_mdlc_random'] == True) | (df['j-is_mdlc_random'] == True)), :].shape[0]

            n_d_prp19 = df.loc[((df['i-is_prp19_random'] == True) | (df['j-is_prp19_random'] == True)), :].shape[0]
            n_d_prp19_metric = df.loc[(df['is_metric'] == True) & ((df['i-is_prp19_random'] == True) | (df['j-is_prp19_random'] == True)), :].shape[0]
            n_d_prp19_ortho = df.loc[(df['is_ortho'] == is_metric_ortho_string) & ((df['i-is_prp19_random'] == True) | (df['j-is_prp19_random'] == True)), :].shape[0]

            r.append((run, n_d_mdlc, n_d_mdlc_metric, n_d_mdlc_ortho, n_d_prp19, n_d_prp19_metric, n_d_prp19_ortho))

        if layer == 'HS':
            random_rnf113_gene_ids = set(random.sample(gene_ids, k=n_is_rnf113_dge))

            df['i-is_rnf113_random'] = df['i'].map(lambda x: x in random_rnf113_gene_ids)
            df['j-is_rnf113_random'] = df['j'].map(lambda x: x in random_rnf113_gene_ids)

            n_d_rnf113 = df.loc[((df['i-is_rnf113_random'] == True) | (df['j-is_rnf113_random'] == True)), :].shape[0]
            n_d_rnf113_metric = df.loc[(df['is_metric'] == True) & ((df['i-is_rnf113_random'] == True) | (df['j-is_rnf113_random'] == True)), :].shape[0]
            n_d_rnf113_ortho = df.loc[(df['is_ortho'] == is_metric_ortho_string) & ((df['i-is_rnf113_random'] == True) | (df['j-is_rnf113_random'] == True)), :].shape[0]

            r.append((run, n_d_rnf113, n_d_rnf113_metric, n_d_rnf113_ortho))

    if layer == 'DM':
        columns = ['run', 'n-d-mdlc', 'n-d-mdlc-metric', 'n-d-mdlc-ortho', 'n-d-prp19', 'n-d-prp19-metric', 'n-d-prp19-ortho']
    elif layer == 'HS':
        columns = ['run', 'n-d-rnf113', 'n-d-rnf113-metric', 'n-d-rnf113-ortho']
    dfr = pd.DataFrame(r, columns=columns)
    print(dfr.agg(['mean', 'std']))
