# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and calculates page rank.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer


if __name__ == '__main__':

    #celltype = args.celltype  # spermatocyte or enterocyte
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #layer = args.layer
    
    for celltype in ['spermatocyte', 'enterocyte']:

        print('Reading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
        rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        G = nx.read_gpickle(rGfile_gpickle)

        for layer in ['HS', 'MM', 'DM']:

            print('Isolate {layer:s} Layer'.format(layer=layer))
            Gt = get_network_layer(G, layer=layer)

            print('Calculating: page rank')
            dict_page_rank = nx.pagerank(Gt, weight='weight')
            dict_page_rank['FBgn0261599']

            dfR = pd.DataFrame({'page_rank': dict_page_rank})

            ##
            # Export
            ##
            print('Saving results to .CSV')
            wCSVFile = 'results/pagerank/{celltype:s}/pagerank-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            ensurePathExists(wCSVFile)
            dfR.to_csv(wCSVFile)

    print('Done.')