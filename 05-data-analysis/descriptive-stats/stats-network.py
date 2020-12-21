# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and prints information.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, get_network_by_attribute, ensurePathExists
from tabulate import tabulate
from itertools import combinations


def df2md(df, y_index=False, *args, **kwargs):
    blob = tabulate(df, headers='keys', tablefmt='pipe', *args, **kwargs)
    if not y_index:
        return '\n'.join(['| {}'.format(row.split('|', 2)[-1]) for row in blob.split('\n')])
    return blob


if __name__ == '__main__':

    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    #
    # Node/Edge stats on complete network
    #
    r = []
    for celltype in celltypes:

        print('Loading {celltype:s} network'.format(celltype=celltype))
        rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}.gpickle'.format(celltype=celltype, network='full')
        G = nx.read_gpickle(rGfile_gpickle)

        for layer in ['HS', 'MM', 'DM']:
            print('Separate {layer:s} layer'.format(layer=layer))
            Gt = get_network_layer(G, layer)

            # Number of nodes/edges
            n_nodes = Gt.number_of_nodes()
            n_edges = Gt.number_of_edges()

            r.append((celltype, layer, n_nodes, n_edges))

    print('# Number of nodes/edges in each layer of the full network\n')
    df_stat = pd.DataFrame(r, columns=['celltype', 'species', '#-nodes', '#-edges'])
    print(df2md(df_stat, floatfmt='.4f'))
    file = 'results/stats-full-network.csv'
    ensurePathExists(file)
    df_stat.to_csv(file)

    #
    # Node/Edge stats on threshold/conserved Network
    #
    network = 'conserved'  # ['thr', 'conserved']
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    if network == 'conserved':
        celltypes = ['spermatocyte', 'enterocyte']
    else:
        celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    #
    r = []
    for celltype in celltypes:

        print('Reading {celltype:s}-{network:s}-{threshold:s} network'.format(celltype=celltype, network=network, threshold=threshold_str))
        path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
        if network in ['thr', 'conserved']:
            rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        elif network == 'full':
            rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}.gpickle'.format(celltype=celltype, network=network)
        G = nx.read_gpickle(rGfile_gpickle)

        for layer in ['HS', 'MM', 'DM']:
            print('Separate {layer:s} layer'.format(layer=layer))
            Gt = get_network_layer(G, layer)

            # Number of nodes/edges
            n_nodes = Gt.number_of_nodes()
            n_edges = Gt.number_of_edges()
            # Number of components (islands)
            n_components = nx.number_connected_components(Gt)
            # Largest Component
            Gtlc = max(nx.connected_component_subgraphs(Gt), key=len)
            n_nodes_largest_component = Gtlc.number_of_nodes()
            n_edges_largest_component = Gtlc.number_of_edges()

            for weight in ['combined_score', 'textmining', 'database', 'experiments', 'coexpression', 'cooccurence', 'fusion', 'neighborhood']:

                edges = [(i, j) for i, j, d in Gt.edges(data=True) if weight in d]
                Gtw = Gt.edge_subgraph(edges).copy()
                # Number of nodes/edges
                n_nodes = Gtw.number_of_nodes()
                n_edges = Gtw.number_of_edges()
                # Number of components (islands)
                n_components = nx.number_connected_components(Gtw)
                # Largest Component
                if n_edges > 0:
                    Gtlc = max(nx.connected_component_subgraphs(Gtw), key=len)
                    n_nodes_largest_component = Gtlc.number_of_nodes()
                    n_edges_largest_component = Gtlc.number_of_edges()
                else:
                    n_nodes_largest_component = 0
                    n_edges_largest_component = 0

                r.append((celltype, layer, weight, n_nodes, n_edges, n_components, n_nodes_largest_component, n_edges_largest_component))

    print('# Number of nodes/edges in the layer of the thresholded>0.5 network\n')

    df_stat = pd.DataFrame(r, columns=['celltype', 'species', 'edge-type', '#-nodes', '#-edges', '#-comps.', '#-nodes-in-lgt-comp.', '#-edges-lgt-comp.'])
    print(df2md(df_stat, floatfmt='.4f'))
    if network in ['thr', 'conserved']:
        file = 'results/stats-{network:s}-{threshold:s}-network.csv'.format(network=network, threshold=threshold_str)
    elif network == 'full':
        file = 'results/stats-{network:s}-network.csv'.format(network=network)
    ensurePathExists(file)
    df_stat.to_csv(file)

    #
    # Pairwise conserved genes
    #
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    biotype = 'protein_coding'

    r = []
    for celltype in celltypes:
        #
        print('Loading {celltype:s}-{network:s} Network'.format(celltype=celltype, network=network))
        path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
        if network == 'thr':
            rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        elif network == 'full':
            rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}.gpickle'.format(celltype=celltype, network=network)
        G = nx.read_gpickle(rGfile_gpickle)

        print('Separate Layers')
        HSG = get_network_layer(G, 'HS')
        MMG = get_network_layer(G, 'MM')
        DMG = get_network_layer(G, 'DM')

        print("Select nodes where biotype='{biotype:s}'".format(biotype=biotype))
        HSG = get_network_by_attribute(HSG, attribute='biotype', value=biotype)
        MMG = get_network_by_attribute(MMG, attribute='biotype', value=biotype)
        DMG = get_network_by_attribute(DMG, attribute='biotype', value=biotype)

        # Pairs
        for (layer_i, Gi), (layer_j, Gj) in combinations([('HS', HSG), ('MM', MMG), ('DM', DMG)], 2):
            print("Comparing: {layer_i:s} with {layer_j:s}".format(layer_i=layer_i, layer_j=layer_j))
            pair = layer_i + 'x' + layer_j

            genes_i = [*Gi.nodes()]
            genes_j = [*Gj.nodes()]

            genes_ij = genes_i + genes_j

            # Only genes between these species
            Gtmp = nx.subgraph(G, genes_ij).copy()

            # Remove intra edges
            remove_intra_edges = [(i, j) for i, j, d in Gtmp.edges(data=True) if d.get('type', None) == 'intra']
            Gtmp.remove_edges_from(remove_intra_edges)

            # Remove isolates
            remove_isolates_nodes = list(nx.isolates(Gtmp))
            Gtmp.remove_nodes_from(remove_isolates_nodes)

            # Keep only homolgos
            Gitmp = nx.subgraph(Gi, Gtmp).copy()
            Gjtmp = nx.subgraph(Gj, Gtmp).copy()

            # Number of nodes/edges
            n_nodes_i = Gitmp.number_of_nodes()
            n_nodes_j = Gjtmp.number_of_nodes()

            n_edges_i = Gitmp.number_of_edges()
            n_edges_j = Gjtmp.number_of_edges()

            r.append((celltype, pair, layer_i, n_nodes_i, n_edges_i, layer_j, n_nodes_j, n_edges_j))

    print('# Pairwise number of conserved nodes/edges of the full network\n')
    df_stat = pd.DataFrame(r, columns=['celltype', 'specie-pair', 'layer-i', '#-nodes-i', '#-edges-i', 'layer-j', '#-nodes-j', '#-edges-j'])
    print(df2md(df_stat, floatfmt='.4f'))
    if network == 'thr':
        file = 'results/stats-{network:s}-{threshold:s}-network.csv'.format(network=network, threshold=threshold_str)
    elif network == 'full':
        file = 'results/stats-{network:s}-pairwise-conserved-network.csv'.format(network=network)
    ensurePathExists(file)
    df_stat.to_csv(file)

