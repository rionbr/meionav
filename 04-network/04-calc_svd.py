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
from scipy import linalg
import networkx as nx
from utils import ensurePathExists, get_network_layer
# from sklearn.preprocessing import Normalizer, PowerTransformer
# from sklearn.decomposition import TruncatedSVD
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
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
    rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    # SVD per Layer
    for layer in ['HS', 'MM', 'DM']:
        print('Isolate {layer:s} Layer'.format(layer=layer))
        Gt = get_network_layer(G, layer=layer)
        #
        dfG = pd.DataFrame(data={'gene': [d.get('label', None) for n, d in Gt.nodes(data=True)]}, index=Gt.nodes)
        #
        if layer == 'DM':
            core = nx.get_node_attributes(Gt, name='core')
            dfG['core'] = dfG.index.map(core)
            #
            fert = nx.get_node_attributes(Gt, name='mean-fert-rate')
            dfG['mean-fert-rate'] = dfG.index.map(fert)
        elif layer in ['MM', 'HS']:
            mammals = nx.get_node_attributes(Gt, name='mammals')
            dfG['mammals'] = dfG.index.map(mammals)

        print('Extract Adjacency Matrix')
        M = nx.to_numpy_matrix(Gt)

        #print('Normalizing')
        #M = Normalizer(norm='l2').fit_transform(M)

        """
        print('Calculating SVD (sklearn)')
        tsvd = TruncatedSVD(n_components=components, algorithm='arpack')
        res = tsvd.fit_transform(M)
        #
        columns = ['{:d}c'.format(i) for i in range(1, components + 1)]
        df_SVD = pd.DataFrame(res[:, 0:components], columns=columns, index=dfG.index)
        df_SVD = pd.concat([dfG, df_SVD], axis='columns')
        #
        s_Var = pd.Series(tsvd.explained_variance_ratio_, index=range(1, (res.shape[1] + 1)), name='explained_variance_ratio')
        """

        print('Calculating SVD (scipy)')
        U, Sigma, V = linalg.svd(M, full_matrices=False, overwrite_a=False, check_finite=False, lapack_driver='gesdd')
        res = U * Sigma
        full_var = np.var(M, axis=0).sum()
        expl_var = np.var(res, axis=0)
        expl_var_ratio = expl_var / full_var

        columns = ['{:d}c'.format(i) for i in range(1, components + 1)]
        df_SVD = pd.DataFrame(res[:, 0:components], columns=columns, index=dfG.index)
        df_SVD = pd.concat([dfG, df_SVD], axis='columns')
        s_Var = pd.Series(expl_var_ratio, index=range(1, (expl_var_ratio.shape[0] + 1)), name='explained_variance_ratio')

        print('Saving results to .CSV')
        wSVDFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        wSFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-s.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        #
        ensurePathExists(wSVDFile)
        ensurePathExists(wSFile)
        #
        df_SVD.to_csv(wSVDFile)
        s_Var.to_csv(wSFile, header=True)

    print('Done.')
