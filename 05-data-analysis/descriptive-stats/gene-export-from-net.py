# coding=utf-8
# Author: Rion B Correia
# Date: Aug 20, 2020
#
# Description: Loads graphs and exports the list of genes
#
# Instructions:
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, ensurePathExists


if __name__ == '__main__':

    celltypes = ['spermatocyte', 'enterocyte', 'neuron', 'muscle']
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    layers = ['HS', 'MM', 'DM']

    for celltype in celltypes:
        print('Loading network: {celltype:s} {network:s} {threshold:s}'.format(celltype=celltype, network=network, threshold=threshold_str))

        path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
        rG_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        rGc_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)

        print('Load graph')
        G = nx.read_gpickle(rG_file_gpickle)
        Gc = nx.read_gpickle(rGc_file_gpickle)

        # Add conserved information on G
        conserved = {node: True for node in Gc.nodes()}
        nx.set_node_attributes(G, name='conserved', values=conserved)

        for layer in layers:

            print('Separate layer {layer:s}'.format(layer=layer))
            Gl = get_network_layer(G, layer)

            df = pd.DataFrame.from_dict(dict(Gl.nodes(data=True)), orient='index')

            # Export
            wCSVFile = 'results/gene-net-export-thr/{celltype:s}/csv-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            ensurePathExists(wCSVFile)
            df.to_csv(wCSVFile)
