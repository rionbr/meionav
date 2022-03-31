# coding=utf-8
# Author: Rion B Correia
# Date: March 30, 2021
#
# Description: Calculates networks backbone stats
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from itertools import product, chain
from collections import Counter

if __name__ == '__main__':

    celltype = 'spermatocyte'
    layers = ['HS', 'MM', 'DM']
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
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    path_net = '../../04-network/results/network/{celltype:}/'.format(celltype=celltype)
    rGfile = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile)
    
    is_not_cross_edges = [(i, j) for i, j, d in G.edges(data=True) if d.get('type') != 'cross']
    H = G.copy()
    H.remove_edges_from(is_not_cross_edges)

    #
    # Conserved
    #
    """
    path_fpkm = '../../02-core_genes/results/FPKM/'
    df_HS = pd.read_csv(path_fpkm + 'HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
    df_MM = pd.read_csv(path_fpkm + 'MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
    df_DM = pd.read_csv(path_fpkm + 'DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')

    dict_string_gene_HS = df_HS['id_gene'].to_dict()
    dict_string_gene_MM = df_MM['id_gene'].to_dict()
    dict_string_gene_DM = df_DM['id_gene'].to_dict()

    print('Loading {celltype:s} meta genes'.format(celltype=celltype))
    path = '../../02-core_genes/results/'
    dfM = pd.read_csv(path + 'meta-genes/meta-{celltype:s}-genes.csv.gz'.format(celltype=celltype), index_col='id_eggnog', usecols=['id_eggnog', 'id_string_HS', 'id_string_MM', 'id_string_DM'])

    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    dfM['id_gene_HS'] = dfM['id_string_HS'].apply(lambda x: [dict_string_gene_HS[i] for i in x])
    dfM['id_gene_MM'] = dfM['id_string_MM'].apply(lambda x: [dict_string_gene_MM[i] for i in x])
    dfM['id_gene_DM'] = dfM['id_string_DM'].apply(lambda x: [dict_string_gene_DM[i] for i in x])

    dfM = dfM[['id_gene_HS', 'id_gene_MM', 'id_gene_DM']]
    # Only keep meta genes with homologs in all three species
    # 1-species
    dfM1 = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 1]
    # 2-species
    dfM2 = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 2]
    # 3-species
    dfM3 = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 3]

    dfM3tmp = dfM3.explode('id_gene_HS').explode('id_gene_MM').explode('id_gene_DM').dropna()
    dict_conserved = {gene: 'conserved-HS-MM-DM' for gene in dfM3tmp.stack().tolist()}
    # Add info to Multilayer Graph
    nx.set_node_attributes(G, values=dict_conserved, name='conserved')
    """
    #
    # Backbones
    # 
    for layer in layers:
        path_backbone = '../../04-network/results/network-closure/{celltype:}/'.format(celltype=celltype)
        rBfile = path_backbone + 'net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        B = nx.read_gpickle(rBfile)
        is_metric_edges = [(i, j) for i, j, d in B.edges(data=True) if d.get('is_metric') == True]
        B = B.edge_subgraph(is_metric_edges).copy()
        #
        #nx.set_node_attributes(B, values=dict_conserved, name='conserved')
        #
        data['backbones'][layer] = B

    # "Ortho backbone" (backbone edges present in another layer)
    r = []
    for layer in layers:
        print("- Computing {layer:s} ortho backbone".format(layer=layer))
        B = data['backbones'][layer]

        number_of_metric_edges = B.number_of_edges()

        dict_is_metric_ortho = {}
        # Loop all backbone edges
        for (i, j) in B.edges():

            at_least_one_metric_ortho = False
            is_metric_ortho_string = 'is_metric_ortho'
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

                is_metric_ortho = any([BO.has_edge(ortho_i, ortho_j) for ortho_i, ortho_j in combinations])

                if is_metric_ortho:
                    at_least_one_metric_ortho = True
                    is_metric_ortho_string += '-{other_layer:s}'.format(other_layer=other_layer)

            if at_least_one_metric_ortho:
                dict_is_metric_ortho[(i, j)] = is_metric_ortho_string
            else:
                dict_is_metric_ortho[(i, j)] = False

        number_of_metric_ortho_edges = len(dict_is_metric_ortho)
        counter = dict(Counter(dict_is_metric_ortho.values()))
        r.append((celltype, layer, number_of_metric_edges, number_of_metric_ortho_edges, counter))
        #
        nx.set_edge_attributes(B, values=dict_is_metric_ortho, name='is_metric_ortho')

    dfR = pd.DataFrame(r, columns=['celltype', 'layer', 'n-edges-metric', 'n-edges-metric-ortho', 'counter'])

    # Export Stats
    dfR.to_csv('results/ortho-backbone/net-backbone-ortho-stats.csv')

    # Export selected edges
    bool_etal = {
        'ENSG00000153046',  # 'CDYL'
        'ENSG00000152430',  # 'BOLL'
        'ENSG00000109832',  # 'DDX25'
        'ENSG00000071626',  # 'DAZAP1'
    }
    tex_etal = {
        'ENSG00000136478',  # TEX2
        'ENSG00000081019',  # RSBN1
        'ENSG00000104450',  # SPAG1
        'ENSG00000145375',  # SPATA5
        'ENSG00000006282',  # SPATA20
    }

    # Select backbone
    B = data['backbones']['HS']

    # Translate name to id
    dfs = pd.read_csv('data/selected-genes.csv', header=0)
    dict_label2id = {d.get('label'): n for n, d in B.nodes(data=True)}
    dfs['id_gene-HS'] = dfs['gene-HS'].map(dict_label2id)
    genes_etal = dfs['id_gene-HS'].tolist()

    genes_edges = B.edges(genes_etal)
    r = []
    for i, j in genes_edges:
        di = B.nodes[i]
        dj = B.nodes[j]
        dij = B[i][j]
        #
        gene_i = di['label']
        gene_j = dj['label']
        logFPKM_i = di['logFPKM']
        logFPKM_j = dj['logFPKM']
        weight = dij['weight']
        is_metric = dij['is_metric']
        is_ultrametric = dij['is_ultrametric']
        is_metric_ortho = dij['is_metric_ortho']
        #
        r.append((i, j, gene_i, gene_j, logFPKM_i, logFPKM_j, weight, is_metric, is_ultrametric, is_metric_ortho))
    dft0 = pd.DataFrame(r, columns=['id_gene_i', 'id_gene_j', 'gene_i', 'gene_j', 'logFPKM_i', 'logFPKM_j', 'weight', 'is_metric', 'is_ultrametric', 'is_metric_ortho'])
    dft0.to_csv('results/ortho-backbone-subgraph/net-ortho-backbone-HS-subgraph-GENES-etal.csv')
    #
    bool_edges = B.edges(bool_etal)
    r = []
    for i, j in bool_edges:
        di = B.nodes[i]
        dj = B.nodes[j]
        dij = B[i][j]
        #
        gene_i = di['label']
        gene_j = dj['label']
        logFPKM_i = di['logFPKM']
        logFPKM_j = dj['logFPKM']
        weight = dij['weight']
        is_metric = dij['is_metric']
        is_ultrametric = dij['is_ultrametric']
        is_metric_ortho = dij['is_metric_ortho']
        #
        r.append((i, j, gene_i, gene_j, logFPKM_i, logFPKM_j, weight, is_metric, is_ultrametric, is_metric_ortho))
    dft1 = pd.DataFrame(r, columns=['id_gene_i', 'id_gene_j', 'gene_i', 'gene_j', 'logFPKM_i', 'logFPKM_j', 'weight', 'is_metric', 'is_ultrametric', 'is_metric_ortho'])
    dft1.to_csv('results/ortho-backbone-subgraph/net-ortho-backbone-HS-subgraph-BOOL-etal.csv')
    #
    tex_edges = B.edges(tex_etal)
    r = []
    for i, j in tex_edges:
        di = B.nodes[i]
        dj = B.nodes[j]
        dij = B[i][j]
        #
        gene_i = di['label']
        gene_j = dj['label']
        logFPKM_i = di['logFPKM']
        logFPKM_j = dj['logFPKM']
        weight = dij['weight']
        is_metric = dij['is_metric']
        is_ultrametric = dij['is_ultrametric']
        is_metric_ortho = dij['is_metric_ortho']
        #
        r.append((i, j, gene_i, gene_j, logFPKM_i, logFPKM_j, weight, is_metric, is_ultrametric, is_metric_ortho))
    dft2 = pd.DataFrame(r, columns=['id_gene_i', 'id_gene_j', 'gene_i', 'gene_j', 'logFPKM_i', 'logFPKM_j', 'weight', 'is_metric', 'is_ultrametric', 'is_metric_ortho'])
    dft2.to_csv('results/ortho-backbone-subgraph/net-ortho-backbone-HS-subgraph-TEX-etal.csv')
    #
    # Identify the genes on the orthologous backbone (for later GOEA analysis)
    #
    for layer in layers:
        print('> Ortho nodes in layer: {layer:s}'.format(layer=layer))
        B = data['backbones'][layer]
        is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])

        is_ortho_metric_edges = [(i, j) for i, j, d in B.edges(data=True) if d.get('is_metric_ortho') == is_metric_ortho_string]
        is_ortho_metric_nodes = set(list(chain(*is_ortho_metric_edges)))

        Z = B.subgraph(is_ortho_metric_nodes).copy()
        dfZ = pd.DataFrame.from_dict(dict(Z.nodes(data=True)), orient='index')

        # Export
        dfZ.to_csv('results/ortho-backbone/net-ortho-backbone-{layer:s}.csv.gz'.format(layer=layer))
