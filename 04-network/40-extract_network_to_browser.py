# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and extracts subgraphs based on parameters for the networkbrowser.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from itertools import chain, product
from utils import ensurePathExists
from matplotlib import colors
import random
import json


cmap_meanfertrate = colors.LinearSegmentedColormap.from_list(name='cmap-mean-fert-rate', colors=['#d62728', '#1f77b4'], N=256)


def fert_rate_color(x):
    if pd.isnull(x):
        return '#FFFFFF'  # white
    else:
        return colors.to_hex(cmap_meanfertrate(x))  # color


def get_network_layer(G, layer=''):
    return G.subgraph([n for n, d in G.nodes(data=True) if (d.get('layer') == layer)]).copy()


def map_new_range_of_node_values(G, attr=''):
    # Retrieve attributes
    dict_node_attr = nx.get_node_attributes(G, name=attr)
    # Makes a set of values
    attset = set(dict_node_attr.values())
    # Makes a map of new integer range to the size of the set
    dict_map = {k: v for k, v in (zip(attset, range(len(attset))))}
    # Applies the map to the original values
    dict_new_node_attr = {k: dict_map[v] for k, v in dict_node_attr.items()}
    # Sets the new attribute values
    nx.set_node_attributes(G, name=attr, values=dict_new_node_attr)
    #
    print("New number of '{att:s}' attr: {n:,d}".format(att=attr, n=len(attset)))
    return G

if __name__ == '__main__':

    network = 'complete'  # 'complete', 'meiotic-entry', 'meiotic-exit'

    print('Reading Network')
    rGfile_gpickle = 'results/net_{network:s}_mlayer_backbone.gpickle'.format(network=network)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Selecting nodes
    #
    print('Selecting nodes')
    # Select only nodes from a layer
    select_layer_nodes = [n for n, d in G.nodes(data=True) if (
        #(d.get('layer') == 'DM') and
        ((d.get('mammals') == True) or (d.get('core') == True))
    )]
    G = G.subgraph(select_layer_nodes).copy()
    #
    # Selecting edges
    #
    # % sample
    #sampled_edges = random.sample(SG.edges, k=300)
    #SG = SG.edge_subgraph(sampled_edges).copy()

    # Only backbone/experimental evidence
    select_layer_edges = [(i, j) for i, j, d in G.edges(data=True) if (
        # ((d.get('type') == 'intra') and (d.get('combined_score', -1) >= 500))
        ((d.get('type') == 'intra') and (d.get('metric-backbone', None) == True))
        or
        (d.get('type') == 'cross')
    )]
    G = G.edge_subgraph(select_layer_edges).copy()

    # Remove Isolates
    isolates = list(nx.isolates(G))
    print('Removing {n:d} isolated nodes'.format(n=len(isolates)))
    for i in isolates:
        print('> {id:s}: {gene:s}'.format(id=i, gene=G.nodes[i]['label']))
    G.remove_nodes_from(isolates)

    print('Translating Module Numbers')
    G = map_new_range_of_node_values(G, 'modules-HS-mammals-infomap')
    G = map_new_range_of_node_values(G, 'modules-HS-mammals-louvain')
    G = map_new_range_of_node_values(G, 'modules-MM-mammals-infomap')
    G = map_new_range_of_node_values(G, 'modules-MM-mammals-louvain')
    G = map_new_range_of_node_values(G, 'modules-DM-infomap')
    G = map_new_range_of_node_values(G, 'modules-DM-louvain')
    G = map_new_range_of_node_values(G, 'modules-infomap')

    print('Adding additional information')
    # Node Color (baed on Mean Fert-Rate)
    #attr_color = {n: fert_rate_color(v) for n, v in G.nodes.data('mean-fert-rate')}
    #nx.set_node_attributes(SG, values=attr_color, name='color-fert-rate')

    # Size
    #nx.set_node_attributes(SG, values=25, name='size')

    # Add node information
    for n, d in G.nodes(data=True):
        id_gene = n

    # Add edge information (layer)
    for i, j, d in G.edges(data=True):
        ijtype = d['type']
        ilayer = G.nodes[i]['layer']
        jlayer = G.nodes[j]['layer']
        if ilayer < jlayer:
            ilayer, jlayer = jlayer, ilayer
        if ijtype == 'intra':
            G[i][j]['layer'] = ilayer
        else:
            G[i][j]['layer'] = '{:s}-{:s}'.format(ilayer, jlayer)

    ##
    # Export
    ##
    print('Exporting')

    # json
    print('> json')
    jsondata = {
        'directed': False,
        'graph': [],
        'nodes': [{'id': i, **d} for i, d in G.nodes(data=True)],
        'edges': [{'from': i, 'to':j, **d} for i, j, d in G.edges(data=True)]
    }
    wSGfile_json = 'results/net_{network:s}_core_mlayer.json'.format(network=network)
    with open(wSGfile_json, 'w') as outfile:
        json.dump(jsondata, outfile, indent=4)

    """
    # graphml
    print('> graphml (multi layer)')
    wSGfile_graphml = 'results/net_core_mlayer.graphml'
    ensurePathExists(wSGfile_graphml)
    nx.write_graphml(SG, wSGfile_graphml)
    # graphml
    print('> graphml (HS layer')
    wHSGfile_graphml = 'results/net_core_HS_layer.graphml'
    ensurePathExists(wHSGfile_graphml)
    nx.write_graphml(get_network_layer(SG, 'HS'), wHSGfile_graphml)

    print('> graphml (MM layer)')
    wMMGfile_graphml = 'results/net_core_MM_layer.graphml'
    ensurePathExists(wMMGfile_graphml)
    nx.write_graphml(get_network_layer(SG, 'MM'), wMMGfile_graphml)

    print('> graphml (DM layer)')
    wDMGfile_graphml = 'results/net_core_DM_layer.graphml'
    ensurePathExists(wDMGfile_graphml)
    nx.write_graphml(get_network_layer(SG, 'DM'), wDMGfile_graphml)
    """