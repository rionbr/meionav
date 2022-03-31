
# coding=utf-8
# Author: Rion B Correia
# Date: March 22, 2021
#
# Description: Reads a MultiLayer module subgraph networks and extracts each layer independently for plotting adding information to color the nodes.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, ensurePathExists
from itertools import combinations, product, chain


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

    #
    # Add Backbone
    #
    path_backbone = '../../04-network/results/network-closure/{celltype:}/'.format(celltype=celltype)
    rBfile_gpickle = path_backbone + 'net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    B = nx.read_gpickle(rBfile_gpickle)
    #
    distance = nx.get_edge_attributes(B, 'distance')
    is_metric = nx.get_edge_attributes(B, 'is_metric')
    is_ultrametric = nx.get_edge_attributes(B, 'is_ultrametric')
    nx.set_edge_attributes(G, values=distance, name='distance')
    nx.set_edge_attributes(G, values=is_metric, name='is_metric')
    nx.set_edge_attributes(G, values=is_ultrametric, name='is_ultrametric')

    to_remove = {(i, j) for i, j, d in B.edges(data=True) if d.get('is_metric') is not True}
    B.remove_edges_from(to_remove)
    #
    # Add Ortho-backbone
    #
    is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])

    path_ortho_backbone = '../../04-network/results/network-closure-ortho/{celltype:s}/'.format(celltype=celltype)
    rOfile_gpickle = path_ortho_backbone + 'net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    O = nx.read_gpickle(rOfile_gpickle)
    #
    is_metric_ortho = nx.get_edge_attributes(O, 'is_metric_ortho')
    nx.set_edge_attributes(G, values=is_metric_ortho, name='is_metric_ortho')

    to_remove = {(i, j) for i, j, d in O.edges(data=True) if d.get('is_metric_ortho', None) != is_metric_ortho_string}
    O.remove_edges_from(to_remove)

    # only layer
    Gt = get_network_layer(G, layer)


    # DEBUG : search for a gene id
    """
    for n,d in M.nodes(data=True):
        if d.get('label') == 'Rbx1-ps':
            print(n,d)
            break
    """

    # General Colors
    colors = {
        'red': '#d62728',
        'orange': '#ff7f0e',
        'pink': '#f7b6d2',
        'other': '#e1e0e2',
    }

    #
    # Ubi
    #
    if layer == 'DM':

        list_ubi = ['FBgn0003943']  # Ubi-p63E
        list_cdc = ['FBgn0036516']  # CG7656
        list_skp = ['FBgn0025637']  # SkpA
        list_roc = ['FBgn0025638']  # Roc1a
        list_cul = ['FBgn0015509']  # Cul1
        list_wee = ['FBgn0011737']  # Wee1
        list_genes = [list_ubi, list_cdc, list_skp, list_roc, list_cul, list_wee]
        #
        node_rename = {
            'CG7656': 'dCdc34',
            'Roc1a': 'dRbx1',
            'SkpA': 'dSkp1',
            'fzy': 'dCdc20',
        }
        node_color = {
            # Source
            'Ubi-p63E': 'red',
            # Intermediate
            'dCdc34': 'orange',
            'dRbx1': 'orange',
            'Cul1': 'orange',
            'dSkp1': 'orange',
            # Target
            'Wee1': 'pink'
        }
        highlight_path = [
            'FBgn0003943',  # Ubi-p63E
            'FBgn0036516',  # CG7656
            'FBgn0025638',  # Roc1a
            'FBgn0015509',  # Cul1
            'FBgn0025637',  # SkpA
            'FBgn0037236',  # Skp2
            'FBgn0001086',  # Fzy
            'FBgn0004106',  # Cdk1
            'FBgn0011737',  # Wee1
        ]
        highlight_edges = set([frozenset([i, j]) for i, j in zip(highlight_path[:-1], highlight_path[1:])])

    elif layer == 'MM':
        list_ubi = ['ENSMUSG00000008348', 'ENSMUSG00000019505']  # Ubc, Ubb
        list_cdc = ['ENSMUSG00000020307', 'ENSMUSG00000020870', 'ENSMUSG00000036241']  # Cdc34, Cdc34b, Ube2r2
        list_skp = ['ENSMUSG00000036309']  # Skp1
        list_roc = ['ENSMUSG00000049832']  # Rbx1-ps
        list_cul = ['ENSMUSG00000029686']  # Cul1
        list_wee = ['ENSMUSG00000031016']  # Wee1
        list_genes = [list_ubi, list_cdc, list_skp, list_roc, list_cul, list_wee]
        node_rename = {
            'Skp1a': 'Skp1',
            'Ube2r2': 'Cdc34b',
            'Cdc34b': 'Cdc34-ps',
            #
            'Rbx1-ps': 'Rbx1',
        }
        node_color = {
            # Source
            'Ubb': 'red',
            'Ubc': 'red',
            'Ube2r2': 'red',
            # Intermediate
            'Cdc34': 'orange',
            'Cdc34b': 'orange',
            'Cdc34-ps': 'orange',
            'Rbx1': 'orange',
            'Cul1': 'orange',
            'Skp1': 'orange',
            # Target
            'Wee1': 'pink'
        }
        highlight_path = [
            'ENSMUSG00000019505',  # Ubb
            'ENSMUSG00000068749',  # Psma4
            'ENSMUSG00000026914',  # Psmd14
            'ENSMUSG00000008348',  # Ubc
            'ENSMUSG00000020307',  # Cdc34
            'ENSMUSG00000029686',  # Cul1
            'ENSMUSG00000036309',  # Skp1a
            'ENSMUSG00000054115',  # Skp2
            'ENSMUSG00000006398',  # Cdc20
            'ENSMUSG00000019942',  # Cdk1
            'ENSMUSG00000044201',  # Cdc25c
            'ENSMUSG00000031016',  # Wee1
        ]
        highlight_edges = set([frozenset([i, j]) for i, j in zip(highlight_path[:-1], highlight_path[1:])])

        # Ubc <-> Cdc34-ps (Cdc34b)
        highlight_edges.add(frozenset(['ENSMUSG00000008348', 'ENSMUSG00000020870']))
        # Ubc <-> Cdc34b (Ube2r2)
        highlight_edges.add(frozenset(['ENSMUSG00000008348', 'ENSMUSG00000036241']))

        # Cdc34-ps (Cdc34b) <-> Rbx1
        highlight_edges.add(frozenset(['ENSMUSG00000020870', 'ENSMUSG00000049832']))
        # Cdc34b (Ube2r2) <-> Rbx1
        highlight_edges.add(frozenset(['ENSMUSG00000036241', 'ENSMUSG00000049832']))

        # Rbx1 <-> Cul1
        highlight_edges.add(frozenset(['ENSMUSG00000049832', 'ENSMUSG00000029686']))

        # Rbx1 <-> Skp1a
        highlight_edges.add(frozenset(['ENSMUSG00000049832', 'ENSMUSG00000036309']))

    elif layer == 'HS':
        list_ubi = ['ENSG00000150991', 'ENSG00000170315', 'ENSG00000263563']  # UBC, UBB, UBBP4
        list_cdc = ['ENSG00000099804', 'ENSG00000107341']  # CDC34, UBE2R2
        list_skp = ['ENSG00000113558']  # Skp1
        list_roc = ['ENSG00000100387']  # RBX1
        list_cul = ['ENSG00000055130']  # CUL1
        list_wee = ['ENSG00000166483']  # WEE1
        list_genes = [list_ubi, list_cdc, list_skp, list_roc, list_cul, list_wee]
        node_rename = {
            'UBE2R2': 'CDC34B'
        }
        node_color = {
            # Source
            'UBB': 'red',
            'UBC': 'red',
            'UBBP4': 'red',
            # Intermediate
            'CDC34': 'orange',
            'CDC34B': 'orange',
            'RBX1': 'orange',
            'CUL1': 'orange',
            'SKP1': 'orange',
            # Target
            'WEE1': 'pink',
        }
        highlight_path = [
            'ENSG00000170315',  # UBB
            'ENSG00000150991',  # UBC
            'ENSG00000099804',  # CDC34
            'ENSG00000100387',  # RBX1
            'ENSG00000055130',  # CUL1
            'ENSG00000113558',  # SKP1
            'ENSG00000173207',  # CKS1B
            'ENSG00000170312',  # CDK1
            'ENSG00000166483',  # WEE1
        ]
        highlight_edges = set([frozenset([i, j]) for i, j in zip(highlight_path[:-1], highlight_path[1:])])
        # Additional edges
        # UBC <-> UBE2R2
        highlight_edges.add(frozenset(['ENSG00000150991', 'ENSG00000107341']))
        # UB2R2 <-> RBX1
        highlight_edges.add(frozenset(['ENSG00000107341', 'ENSG00000100387']))
        # CDK1 <-> CDC20
        highlight_edges.add(frozenset(['ENSG00000170312', 'ENSG00000117399']))
        # CDC34 <-> CUL1
        highlight_edges.add(frozenset(['ENSG00000099804', 'ENSG00000055130']))
        # RBX1 <-> SKP1
        highlight_edges.add(frozenset(['ENSG00000100387', 'ENSG00000113558']))

    all_ids = set(chain(*list_genes))
    edges = set()
    #
    for list_a, list_b in combinations(list_genes, 2):

        for a, b in product(list_a, list_b):
            # direct path
            #path = nx.shortest_path(Gt, source=a, target=b, weight=None)
            #for i, j in zip(path[:-1], path[1:]):
            #    edges.add((i, j))

            # backbone path
            path = nx.shortest_path(B, source=a, target=b, weight='metric_distance')
            for i, j in zip(path[:-1], path[1:]):
                edges.add((i, j))

    #Gt['FBgn0004106']['FBgn0011737']
    #[Gt.nodes[n]['label'] for n in path]

    # Compose the two subgraphs
    Gtf = nx.edge_subgraph(Gt, edges).copy()

    # String has an ID issue with "Rbx1-ps" and "Rbx1". In the network, "Rbx1-ps" is actually "Rbx1"
    if layer =='MM':
        rbx1ps = 'ENSMUSG00000049832'
        #rbx1 = 'ENSMUSG00000022400'
        Gtf.nodes[rbx1ps]['label'] = 'Rbx1'
        Gtf.nodes[rbx1ps]['TPM'] = 29.239579000000006
        Gtf.nodes[rbx1ps]['logFPKM'] = 4.216756513893444
        Gtf.nodes[rbx1ps]['biotype'] = 'protein_coding'

    # Rename Nodes
    for i, d in Gtf.nodes(data=True):
        label = d.get('label', None)
        if label in node_rename.keys():
            new_label = node_rename[label]
            Gtf.nodes[i]['label'] = new_label

    # Color nodes
    for i, d in Gtf.nodes(data=True):
        label = d.get('label', None)
        if label in node_color:
            color = colors[node_color[label]]
        else:
            color = colors['other']
        Gtf.nodes[i]['color'] = color

    # Color edges
    for i, j, d in Gtf.edges(data=True):
        if d.get('is_metric_ortho', None) == is_metric_ortho_string:
            if set([i, j]) in highlight_edges:
                Gtf[i][j]['color'] = '#d62728'  # dark red
            else:
                Gtf[i][j]['color'] = '#f6d3d4'  # light red
        elif d.get('is_metric', None) is True:
            if set([i, j]) in highlight_edges:
                Gtf[i][j]['color'] = '#2ca02c'  # dark green
            else:
                Gtf[i][j]['color'] = '#d4ecd4'  # light green
        else:
            Gtf[i][j]['color'] = '#c7c7c7'  # light gray

        Gtf[i][j]['weight'] = Gtf[i][j]['weight'] * 10 / 2


    # Export
    wGefile_graphml = 'results/gene-UBI-2-WEE/net-path-ubi-2-wee-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    ensurePathExists(wGefile_graphml)
    nx.write_graphml(Gtf, wGefile_graphml)
