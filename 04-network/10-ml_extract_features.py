# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and extracts features for a Machine Learning experiment (extracts X).
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, ensurePathExists
from distanceclosure.utils import _prox2dist as prox2dist
from sklearn.decomposition import PCA
from data import *
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    layers = ['HS', 'MM', 'DM']
    data_cell = data_cells[celltype]

    #
    print('Reading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
    rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    for layer in layers:

        print('Extracting {layer:s} SubGraph'.format(layer=layer))
        Gt = get_network_layer(G, layer)
        #data_cell['graphs'][layer] = Gt

        print('Extracting {layer:s} modules'.format(layer=layer))
        modules = data_cell['modules-svd']['modules'][layer]

        print('Loading SVD for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
        rPCAFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        dfPCA = pd.read_csv(rPCAFile, index_col=0)

        print('Loading Backbone for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
        rBNFile = 'results/backbone/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-backbone.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        dfBN = pd.read_csv(rBNFile, index_col=['i', 'j'], header=0, names=['i', 'j', 'backbone', 's_values'], dtype={'backbone': float, 's_values': float})
        dfBN['backbone'] = dfBN['backbone'].astype(bool)
        dfBN = dfBN.loc[dfBN['backbone'] == True, :]
        backbone_edge_ids = dfBN['backbone'].to_dict()
        nx.set_edge_attributes(Gt, name='metric-backbone', values=backbone_edge_ids)

        # Selecting SVD Component
        for module in modules:
            mid = module['id']
            name = module['name']
            
            print("Extracting Features for module {mid:}-{name:}".format(mid=mid, name=name))
            #
            icx = "{:d}c".format(module['xy-coords']['x-comp'])
            icy = "{:d}c".format(module['xy-coords']['y-comp'])
            icxl, icxh = module['xy-coords']['x-values']
            icyl, icyh = module['xy-coords']['y-values']

            dfPCA_tmp = dfPCA.loc[
                (
                    (dfPCA[icx] >= icxl) & (dfPCA[icx] <= icxh) &
                    (dfPCA[icy] >= icyl) & (dfPCA[icy] <= icyh)
                ), ['gene', icx, icy]]

            component_node_ids = dfPCA_tmp.index.to_list()

            #
            Gtc = nx.subgraph(Gt, component_node_ids)

            # Converting prox2dist
            nx.set_edge_attributes(Gtc, name='distance', values={(i, j): prox2dist(d['weight']) for i, j, d in Gt.edges(data=True)})
            #
            # Caculate Network Features
            #
            print('Calculating: degree centrality')
            dict_degree_centrality = nx.degree_centrality(Gt)
            nx.set_node_attributes(Gtc, name='degree_centrality', values=dict_degree_centrality)

            print("Calculating: eigevector centrality")
            dict_eigenvector_centrality = nx.eigenvector_centrality(Gtc, weight='weight')
            nx.set_node_attributes(Gtc, name='eigenvector_centrality', values=dict_eigenvector_centrality)

            #print('Calculating: katz centrality')
            #dict_katz_centrality = nx.katz_centrality(Gt, weight='weight')
            #nx.set_node_attributes(Gt, name='katz_centrality', values=dict_katz_centrality)

            print('Calculating: closeness centrality')
            dict_closeness_centrality = nx.closeness_centrality(Gtc, distance='distance')

            print('Calculating: betweenness centrality (this may take a while)')
            dict_betweenness_centrality = nx.betweenness_centrality(Gtc)
            nx.set_node_attributes(Gtc, name='betweenness_centrality', values=dict_betweenness_centrality)

            print('Calculating: page rank')
            dict_page_rank = nx.pagerank(Gtc)
            nx.set_node_attributes(Gtc, name='pagerank', values=dict_page_rank)

            print('Calculating: clustering (this may take a while)')
            dict_clustering = nx.clustering(Gtc, weight='weight')
            nx.set_node_attributes(Gtc, name='clustering', values=dict_clustering)

            print('Calculating: average neighbor degree')
            dict_average_neighbor_degree = nx.average_neighbor_degree(Gtc, weight='weight')
            nx.set_node_attributes(Gtc, name='average_neighbor_degree', values=dict_average_neighbor_degree)

            print('Calculating: k nearest neighbors')
            dict_k_nearest_neighbors = nx.k_nearest_neighbors(Gtc, weight='weight')
            nx.set_node_attributes(Gtc, name='k_nearest_neighbors', values=dict_k_nearest_neighbors)

            print('Calculating: eccentricity')
            dict_eccentricity = nx.eccentricity(Gtc)
            nx.set_node_attributes(Gtc, name='eccentricity', values=dict_eccentricity)

            print('Calculating: core number')
            dict_core_number = nx.core_number(Gtc)
            nx.set_node_attributes(Gtc, name='core_number', values=dict_core_number)

            print('Calculating: degree')
            dict_degree = {n: d for n, d in list(Gtc.degree())}
            nx.set_node_attributes(Gtc, name='degree', values=dict_degree)
            dict_degree_weight = {n: d for n, d in list(Gtc.degree(weight='weight'))}
            nx.set_node_attributes(Gtc, name='degree_weight', values=dict_degree_weight)
            dict_degree_backbone = {n: len([n for i, j, d in Gtc.edges(n, data=True) if d.get('metric-backbone', False) == True]) for n in Gtc.nodes()}
            nx.set_node_attributes(Gtc, name='degree_backbone', values=dict_degree_backbone)
            dict_degree_experiments = {n: d for n, d in list(Gtc.degree(weight='experiments'))}
            nx.set_node_attributes(Gtc, name='degree_experiments_weight', values=dict_degree_experiments)

            print('Calculating: PCA')
            pca = PCA(n_components=9)
            M = nx.to_numpy_matrix(Gtc)
            res = pca.fit_transform(M)
            dfSVD = pd.DataFrame(res[:, 0:9], columns=['PCA-1c', 'PCA-2c', 'PCA-3c', 'PCA-4c', 'PCA-5c', 'PCA-6c', 'PCA-7c', 'PCA-8c', 'PCA-9c'], index=Gtc.nodes())

            # Generate dfX
            dfNet = pd.DataFrame.from_dict(dict(Gtc.nodes(data=True)), orient='index')

            # Merge DataFrames
            dfX = pd.concat([dfNet, dfSVD], axis='columns')

            # Calculate y
            """
            def calc_y(r):
                if r.get('mean-fert-rate', 1.0) < 0.7:
                    return True
                elif not pd.isnull(r.get('known-DM-phenotype', None)):
                    return True
                elif not pd.isnull(r.get('new-DM-phenotype', None)):
                    return True
                else:
                    return False
            dfX['y'] = dfX.apply(calc_y, axis='columns')
            """

            ##
            # Export
            ##
            print('Saving results to .CSV')
            wMLFile = 'results/ml/{celltype:s}/{layer:s}/ml-{celltype:s}-{layer:s}-mod-{mid:d}.csv.gz'.format(celltype=celltype, layer=layer, mid=mid)
            ensurePathExists(wMLFile)
            dfX.to_csv(wMLFile)
