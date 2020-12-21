# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and computed its modules using Louvain & Infomap.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
# from scipy import linalg
# from sklearn.preprocessing import Normalizer, PowerTransformer
# from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import PCA
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    #parser.add_argument("--layer", default='DM', type=str, choices=['DM', 'MM', 'HS'], help="Network layer to compute SVD. Defaults to 'DM'.")
    parser.add_argument("--components", default=15, type=int, help="Number of singular values (components) to calculate. Defaults to 15.")

    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    #layer = args.layer
    components = args.components

    # For "sign indeterminacy"
    np.random.seed(1)
    #
    # Load Network
    #
    print('Reading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
    rGfile_gpickle = 'results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    # SVD per Layer
    for layer in ['HS', 'MM', 'DM']:
        print('Isolate {layer:s} Layer'.format(layer=layer))
        Gt = get_network_layer(G, layer=layer)
        #
        dfG = pd.DataFrame(data={'gene': [d.get('label', None) for n, d in Gt.nodes(data=True)]}, index=Gt.nodes)
        #
        print('Extract Adjacency Matrix')
        M = nx.to_numpy_matrix(Gt)

        print('Calculating PCA (sklearn)')
        pca = PCA(n_components=None, svd_solver='full')
        res = pca.fit(M).transform(M)
        #
        columns = ['{:d}c'.format(i) for i in range(1, components + 1)]
        df_pca = pd.DataFrame(res[:, 0:components], columns=columns, index=dfG.index)
        df_pca = pd.concat([dfG, df_pca], axis='columns')
        #
        s_pca_var = pd.Series(pca.explained_variance_ratio_, index=range(1, (res.shape[1] + 1)), name='explained_variance_ratio')

        print('Saving results to .CSV')
        wPCAFile = 'results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        wSFile = 'results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-s.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        #
        ensurePathExists(wPCAFile)
        ensurePathExists(wSFile)
        #
        df_pca.to_csv(wPCAFile)
        s_pca_var.to_csv(wSFile, header=True)

    print('Done.')
