# coding=utf-8
# Author: Rion B Correia
# Date: Nov 11, 2021
#
# Description: Computes overlaps between DGE results for gene mdlc (DM) and RFG113 (HS).
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

import networkx as nx
from utils import get_network_layer
from itertools import chain


if __name__ == '__main__':

    # DGE
    df_dge_hs = pd.read_csv('../../01-diff-gene-exp/results/rnf113/HS-DGE-rnf113_vs_control-SELECTED.csv', index_col=0)
    df_dge_dm = pd.read_csv('../../01-diff-gene-exp/results/mdlc/DM-DGE-mdlc_vs_control-SELECTED.csv', index_col=0)

    # Core
    df_core_hs = pd.read_csv('../../02-core_genes/results/pipeline-core/HS_meiotic_genes.csv', index_col=0)
    df_core_dm = pd.read_csv('../../02-core_genes/results/pipeline-core/DM_meiotic_genes.csv', index_col=0)
    dict_core_hs = {i: True for i in df_core_hs.index.tolist()}
    dict_core_dm = {i: True for i in df_core_dm.index.tolist()}

    dge_hs = set(df_dge_hs.index.to_list())
    dge_dm = set(df_dge_dm.index.to_list())
    n_dge_hs = len(dge_hs)
    n_dge_dm = len(dge_dm)
    print("n DGE HS: {n_dge_hs:d}".format(n_dge_hs=n_dge_hs))
    print("n DGE DM: {n_dge_dm:d}".format(n_dge_dm=n_dge_dm))

    G = nx.read_gpickle('results/net-spermatocyte-rnf113.gpickle')

    # Thresholded Network
    gfile = '../../04-network/results/network/spermatocyte/net-spermatocyte-thr-0p5.gpickle'
    TG = nx.read_gpickle(gfile)

    layers = ['HS', 'MM', 'DM']
    for layer in ['HS', 'DM']:

        # Layer of the thresholded network
        TGt = get_network_layer(TG, layer=layer)

        # Page Rank
        dict_page_rank = nx.pagerank(TGt)
        nx.set_node_attributes(G, name='page-rank', values=dict_page_rank)
        #

        gfile = '../../04-network/results/network-closure-ortho/spermatocyte/net-closure-ortho-spermatocyte-thr-0p5-{layer:s}.gpickle'.format(layer=layer)
        OB = nx.read_gpickle(gfile)

        #
        is_metric = nx.get_edge_attributes(OB, name='is_metric')
        is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])
        is_metric_ortho = {k: v for k, v in nx.get_edge_attributes(OB, name='is_metric_ortho').items() if v == is_metric_ortho_string}
        # Nodes
        is_ortho_metric_edges = [(i, j) for i, j, d in OB.edges(data=True) if d.get('is_metric_ortho') == is_metric_ortho_string]
        is_ortho_metric_nodes = {i: True for i in chain(*is_ortho_metric_edges)}
        nx.set_node_attributes(G, name='is_metric_ortho', values=is_ortho_metric_nodes)
        # Edges
        nx.set_edge_attributes(G, name='is_metric', values=is_metric)
        nx.set_edge_attributes(G, name='is_metric_ortho', values=is_metric_ortho)

    # Add Module information
    for layer in ['HS', 'DM']:
        dfile = '../entropy-based-modules/results/pca-entropy/spermatocyte/{layer:s}/pca-spermatocyte-thr-0p5-{layer:s}-modules.csv.gz'.format(layer=layer)
        dfm = pd.read_csv(dfile, index_col=0)
        dict_module = dfm['module-id'].to_dict()
        #
        nx.set_node_attributes(G, name='module', values=dict_module)

    # Add Core Attribute
    nx.set_node_attributes(G, name='core', values=dict_core_hs)
    nx.set_node_attributes(G, name='core', values=dict_core_dm)

    n_nodes = G.number_of_nodes()
    n_nodes_hs = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'HS'])
    n_nodes_dm = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'DM'])

    n_nodes_core_hs = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'HS' and d.get('core', None) is True])
    n_nodes_core_dm = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'DM' and d.get('core', None) is True])

    n_nodes_ortho_hs = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'HS' and d.get('is_metric_ortho', None) is True])
    n_nodes_ortho_mm = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'MM' and d.get('is_metric_ortho', None) is True])

    n_nodes_down_hs = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'HS' and d.get('logFC', None) < 0])
    n_nodes_down_dm = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'DM' and d.get('logFC', None) < 0])

    n_nodes_coredown_hs = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'HS' and d.get('core', None) is True and d.get('logFC', None) < 0])
    n_nodes_coredown_dm = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'DM' and d.get('core', None) is True and d.get('logFC', None) < 0])

    n_nodes_orthodown_hs = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'HS' and d.get('logFC', None) < 0 and d.get('is_metric_ortho', None) is True])
    n_nodes_orthodown_dm = len([i for i, d in G.nodes(data=True) if d.get('layer', None) == 'DM' and d.get('logFC', None) < 0 and d.get('is_metric_ortho', None) is True])

    print(">Multilayer mdlc/RNF113 network")

    print("n nodes: {n_nodes:d}".format(n_nodes=n_nodes))

    print("n nodes HS: {n_nodes_hs:d}".format(n_nodes_hs=n_nodes_hs))
    print("n nodes DM: {n_nodes_dm:d}".format(n_nodes_dm=n_nodes_dm))

    print("n nodes core HS: {n_nodes_core_hs:d}".format(n_nodes_core_hs=n_nodes_core_hs))
    print("n nodes core DM: {n_nodes_core_dm:d}".format(n_nodes_core_dm=n_nodes_core_dm))

    print("n nodes down HS: {n_nodes_down_hs:d}".format(n_nodes_down_hs=n_nodes_down_hs))
    print("n nodes down DM: {n_nodes_down_dm:d}".format(n_nodes_down_dm=n_nodes_down_dm))

    print("n nodes core-down HS: {n_nodes_coredown_hs:d}".format(n_nodes_coredown_hs=n_nodes_coredown_hs))
    print("n nodes core-down DM: {n_nodes_coredown_dm:d}".format(n_nodes_coredown_dm=n_nodes_coredown_dm))

    print("n nodes ortho-down HS: {n_nodes_orthodown_hs:d}".format(n_nodes_orthodown_hs=n_nodes_orthodown_hs))
    print("n nodes ortho-down DM: {n_nodes_orthodown_dm:d}".format(n_nodes_orthodown_dm=n_nodes_orthodown_dm))

    dfR = pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient='index').sort_values('layer')

    # how to find all genes simultaneosly downregulated in DM and HS.
    # In case where 1-to-N relations, consider at least one neighbor downreg.
    dict_nei_reg = {}
    for i, di in G.nodes(data=True):
        jlist = []
        for j in nx.neighbors(G, i):
            if G.nodes[i]['layer'] != G.nodes[j]['layer']:
                if G.nodes[j].get('logFC') < 0:
                    jlist.append('down')
                else:
                    jlist.append('up')
        dict_nei_reg[i] = ','.join(jlist)
    dfR['reg'] = dfR['logFC'].map(lambda x: 'up' if x > 0 else 'down')
    dfR['cross-reg'] = dfR.index.map(dict_nei_reg)

    # Add cross neighbor information
    dict_nei = {}
    for i in G.nodes():
        jlist = []
        for j in nx.neighbors(G, i):
            if G.nodes[i]['layer'] != G.nodes[j]['layer']:
                jlist.append(j)
        if len(jlist):
            dict_nei[i] = ','.join(jlist)
    #
    dfR['id-cross'] = dfR.index.map(dict_nei)

    # Select only HS - orthoBB - downreg genes in both species
    dfZ = dfR.loc[(dfR['layer'] == 'HS') & (dfR['is_metric_ortho'] == True) & (dfR['reg'] == 'down') & (dfR['cross-reg'].str.contains('down')), :]

    # Remove some nodes? (Paulo defined as: para ficarem só aqueles cujo o ortólogo está downregulated nas 2 espécies)
    nodes2remove = {
        'ENSG00000166233',
        'ENSG00000120910',
        'ENSG00000222038',
    }
    dfZ = dfZ.loc[~dfZ.index.isin(nodes2remove), :]

    # Subgraph
    Z = nx.subgraph(G, nbunch=dfZ.index).copy()

    # Remove isolates
    Z.remove_nodes_from(list(nx.isolates(Z)))

    # nodes colors
    colors = {
        'highlight': '#f7032d',  # red
        'evidence': '#f96d6b'  # pink
    }
    nhl = {
        # post-meiotic
        'ENSG00000174953': {'name': 'DHX36', 'type': 'evidence'},
        'ENSG00000113300': {'name': 'CNOT6', 'type': 'evidence'},
        'ENSG00000003987': {'name': 'MTMR7', 'type': 'evidence'},
        # meiotic
        'ENSG00000113558': {'name': 'SKP1', 'type': 'evidence'},
        'ENSG00000125818': {'name': 'PSMF1', 'type': 'evidence'},
        'ENSG00000126803': {'name': 'HSPA2', 'type': 'highlight'},
        'ENSG00000078140': {'name': 'UBE2K', 'type': 'evidence'},
        'ENSG00000152430': {'name': 'BOLL', 'type': 'evidence'},
        'ENSG00000079335': {'name': 'CDC14A', 'type': 'evidence'},
        # post-meotic
        'ENSG00000036257': {'name': 'CUL3', 'type': 'evidence'},
        'ENSG00000144895': {'name': 'EIF2A', 'type': 'evidence'},
        'ENSG00000141232': {'name': 'TOB1', 'type': 'evidence'},
        'ENSG00000198791': {'name': 'CNOT7', 'type': 'evidence'},
        'ENSG00000079435': {'name': 'LIPE', 'type': 'evidence'},
        'ENSG00000198604': {'name': 'BAZ1A', 'type': 'evidence'},
        'ENSG00000182481': {'name': 'KPNA2', 'type': 'highlight'},
        # gametes
        'ENSG00000122735': {'name': 'DNAI1', 'type': 'evidence'},
        'ENSG00000171595': {'name': 'DNAI2', 'type': 'evidence'},
        'ENSG00000184203': {'name': 'PPP1R2', 'type': 'evidence'},
        'ENSG00000118997': {'name': 'DNAH7', 'type': 'evidence'},
        'ENSG00000089101': {'name': 'CFAP61', 'type': 'evidence'},
        'ENSG00000127824': {'name': 'TUBA4A', 'type': 'evidence'},
        'ENSG00000124721': {'name': 'DNAH8', 'type': 'evidence'},
        'ENSG00000096063': {'name': 'SRPK1', 'type': 'evidence'},
        'ENSG00000159625': {'name': 'DRC7', 'type': 'evidence'},
        'ENSG00000183833': {'name': 'CFAP91', 'type': 'evidence'},
        # unspecified
        'ENSG00000134779': {'name': 'TPGS2', 'type': 'evidence'},
    }
    # Nodes
    for i, d in Z.nodes(data=True):
        if i in nhl.keys():
            Z.nodes[i]['color'] = colors[nhl[i]['type']]
        else:
            Z.nodes[i]['color'] = '#d9d9d9'  # light gray
    # Edges
    for i, j, d in Z.edges(data=True):
        # EXPERIENCE coloring
        if d.get('experiments', None) is not None:
            Z[i][j]['color'] = '#d692e7'  # pink
        elif d.get('coexpression', None) is not None:
            Z[i][j]['color'] = '#f0bb57'  # orange
        elif d.get('database', None) is not None:
            Z[i][j]['color'] = '#aec7e8'  # blue
        elif d.get('textmining', None) is not None:
            Z[i][j]['color'] = '#a3d6d2'  # cyan
        else:
            Z[i][j]['color'] = '#c7c7c7'


    nx.write_graphml(Z, 'results/net-multilayer-rnf113-hs.graphml')
    # Export CSV
    dfR.to_csv('results/net-multilayer-rnf113.csv')


