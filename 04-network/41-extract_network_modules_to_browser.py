
# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and extracts subgraphs based on parameters for the networkbrowser.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists
from matplotlib import colors
import json
from utils import get_network_layer

cmap_meanfertrate = colors.LinearSegmentedColormap.from_list(name='cmap-mean-fert-rate', colors=['#d62728', '#1f77b4'], N=256)


def fert_rate_color(x):
    if pd.isnull(x):
        return '#FFFFFF'  # white
    else:
        return colors.to_hex(cmap_meanfertrate(x))  # color


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

    network = 'threshold'  # 'complete', 'backbone', 'threshold'
    level = 0.5
    levelstr = str(level).replace('.', 'p')

    print('Reading Network')
    rGfile_gpickle = 'results/net_{network:s}-{level:s}_mlayer.gpickle'.format(network=network, level=levelstr)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Selecting components
    print('Selecting components')

    # this data comes from '04x-set_svd_modules.py'
    data = {
        'HS': [
            {'name': 'Ubiquitination', 'component': '1c'},
            {'name': 'Translation + Splicing', 'component': '2c'},
            {'name': 'Splicing', 'component': '3c'},
            {'name': 'Cell cycle', 'component': '4c'},
            {'name': 'Unnamed Comp.5', 'component': '5c'},
        ],
        'MM': [
            {'name': 'Ubiquitination', 'component': '1c'},
            {'name': 'Translation', 'component': '2c'},
            {'name': 'Splicing', 'component': '3c'},
            {'name': 'Cell Cycle', 'component': '4c'},
            {'name': 'Unnamed Comp.5', 'component': '5c'},
        ],
        'DM': [
            {'name': 'Translation', 'component': '1c'},
            {'name': 'Splicing', 'component': '2c'},
            {'name': 'Ubiquitination + Cell cycle', 'component': '3c'},
            {'name': 'Post-meiotic Development', 'component': '4c'},
            {'name': 'Unnamed Comp.5 > 6', 'component': '5c'},
            {'name': 'Unnamed Comp.6 > 5', 'component': '6c'},
        ]
    }

    for layer, modules in data.items():

        print('Isolate {layer:s} Layer'.format(layer=layer))
        Gt = get_network_layer(G, layer=layer)

        # Extract Component Modules
        for module in modules:

            name = module['name']
            component = module['component']
            attribute = 'modules-{layer:s}-SVD-{component:s}'.format(layer=layer, component=component)

            print('Extracting module: {module:s}'.format(module=name))

            # Select only nodes from the module
            select_component_nodes = [n for n, d in Gt.nodes(data=True) if (d.get(attribute, False) == True)]
            Gtc = Gt.subgraph(select_component_nodes).copy()

            # Remove Isolates
            """
            isolates = list(nx.isolates(Gtc))
            print('Removing {n:d} isolated nodes'.format(n=len(isolates)))
            for i in isolates:
                print('> {id:s}: {gene:s}'.format(id=i, gene=Gtc.nodes[i]['label']))
            Gtc.remove_nodes_from(isolates)
            """

            # Node Color (baed on Mean Fert-Rate)
            print('Adding additional information')
            """
            attr_color = {n: fert_rate_color(v) for n, v in Gtc.nodes.data('mean-fert-rate')}
            nx.set_node_attributes(SG, values=attr_color, name='color-fert-rate')

            # Size
            nx.set_node_attributes(SG, values=25, name='size')
            """

            # Add node information
            for n, d in Gtc.nodes(data=True):

                id_gene = n

            # Add edge information (layer)
            for i, j, d in Gtc.edges(data=True):
                ijtype = d['type']
                ilayer = Gtc.nodes[i]['layer']
                jlayer = Gtc.nodes[j]['layer']
                if ilayer < jlayer:
                    ilayer, jlayer = jlayer, ilayer
                if ijtype == 'intra':
                    Gtc[i][j]['layer'] = ilayer
                else:
                    Gtc[i][j]['layer'] = '{:s}-{:s}'.format(ilayer, jlayer)

            ##
            # Export
            ##
            print('Exporting')

            # json
            print('> json')
            jsondata = {
                'directed': False,
                'graph': [],
                'nodes': [{'id': i, **d} for i, d in Gtc.nodes(data=True)],
                'edges': [{'from': i, 'to': j, **d} for i, j, d in Gtc.edges(data=True)]
            }
            wGtcfile_json = 'results/json/net_{network:s}-{level:s}-{layer:s}-SVD-{component:s}.json'.format(network=network, level=levelstr, layer=layer, component=component)
            ensurePathExists(wGtcfile_json)
            with open(wGtcfile_json, 'w') as outfile:
                json.dump(jsondata, outfile, indent=4)

            """
            # graphml
            print('> graphml')
            wGtcfile_graphml = 'results/graphml/net_{network:s}-{layer:s}-SVD-{component:s}.graphml'.format(network=network, layer=layerstr, component=component)
            ensurePathExists(wGtcfile_graphml)
            nx.write_graphml(Gtc, wDMGfile_graphml)
            """