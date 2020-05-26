# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and computed its backbone.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = 'thr'
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')

    #
    #
    #
    print('Reading Network')
    rGfile_gpickle = 'results/net_{celltype:s}.gpickle'.format(celltype=celltype)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Thresholding
    #
    print('Removing edges (remove edges < threshold)')
    print('Number of edges: {:,d}'.format(G.number_of_edges()))
    edges_to_remove = [(i, j) for i, j, d in G.edges(data=True) if (
        ((d.get('type') == 'intra') and (d.get('weight') < threshold))
    )]
    G.remove_edges_from(edges_to_remove)
    print('Number of edges: {:,d}'.format(G.number_of_edges()))

    ##
    # Export
    ##
    print('Exporting')
    wGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    ensurePathExists(wGfile_gpickle)
    nx.write_gpickle(G, wGfile_gpickle)
