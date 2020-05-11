
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
from utils import get_network_layer
import argparse
from data import *


cmap_meanfertrate = colors.LinearSegmentedColormap.from_list(name='cmap-mean-fert-rate', colors=['#d62728', '#1f77b4'], N=256)


def fert_rate_color(x):
    if pd.isnull(x):
        return '#FFFFFF'  # white
    else:
        return colors.to_hex(cmap_meanfertrate(x))  # color


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    # parser.add_argument("--layer", default='DM', type=str, choices=['DM', 'MM', 'HS'], help="Network layer to compute SVD. Defaults to 'DM'.")

    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    #layer = args.layer

    #
    data_cell = data_cells[celltype]
    #
    print('Reading Network')
    rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    for layer in ['HS', 'MM', 'DM']:

        print('Isolate {layer:s} Layer'.format(layer=layer))
        Gt = get_network_layer(G, layer=layer)

        modules = data_cell['modules-svd']['modules'][layer]

        # Add Module Information
        print('Load SVD Results ({layer:s})'.format(layer=layer))
        #
        rPCAFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        dfPCA = pd.read_csv(rPCAFile, index_col=0, encoding='utf-8')

        # Extract Component Modules
        for module in modules:
            mid = module['id']
            name = module['name']

            cx = "{:d}c".format(module['xy-coords']['x-comp'])
            cy = "{:d}c".format(module['xy-coords']['y-comp'])
            cxl, cxh = module['xy-coords']['x-values']
            cyl, cyh = module['xy-coords']['y-values']
            query = ((dfPCA[cx] >= cxl) & (dfPCA[cx] <= cxh) & (dfPCA[cy] >= cyl) & (dfPCA[cy] <= cyh))
            component_ids = {g: True for g in dfPCA.loc[query, ['gene', cx, cy]].index.tolist()}
            net_attribute_name = 'module-svd-{layer:s}-{mid:d}'.format(layer=layer, mid=mid)
            nx.set_node_attributes(Gt, values=component_ids, name=net_attribute_name)


        # Remove Isolates
        """
        isolates = list(nx.isolates(Gt))
        print('Removing {n:d} isolated nodes'.format(n=len(isolates)))
        Gt.remove_nodes_from(isolates)
        """

        # Largest Connected Component
        Gt = max(nx.connected_component_subgraphs(Gt), key=len)

        # graphml
        print('Export to graphml')
        wGtfile_graphml = 'results/graphml/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        ensurePathExists(wGtfile_graphml)
        nx.write_graphml(Gt, wGtfile_graphml)
