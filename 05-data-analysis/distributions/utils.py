import os
import networkx as nx


def get_network_layer(G, layer=''):
    return G.subgraph([n for n, d in G.nodes(data=True) if (d.get('layer') == layer)]).copy()


def get_network_by_attribute(G, attribute='', value=''):
    return G.subgraph([n for n, d in G.nodes(data=True) if (d.get(attribute) == value)]).copy()


def get_network_largest_connected_component(G):
    largest_cc = max(nx.connected_components(G), key=len)
    return G.subgraph(largest_cc).copy()


def ensurePathExists(path):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
