# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and the SVD results and calculats a similrity between SVD Modules across layers using gene homology.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from itertools import combinations
from utils import ensurePathExists, get_network_layer
from data import *


def calc_module_proximity(celltype='spermatocyte', network='thr', threshold=0.5):
    threshold_str = str(threshold).replace('.', 'p')
    #
    layers = ['HS', 'MM', 'DM']
    data_cell = data_cells[celltype]

    print('Reading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
    rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    for layer in layers:
        print('Isolate {layer:s} Layer'.format(layer=layer))
        Gt = get_network_layer(G, layer=layer)
        data_cell['graphs'][layer] = Gt

        print('Loading SVD for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
        rPCAFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        dfPCA = pd.read_csv(rPCAFile, index_col=0, encoding='utf-8')
        data_cell['modules-svd']['files'][layer] = rPCAFile
        data_cell['modules-svd']['dfs'][layer] = dfPCA

    r = []
    for (layer_i), (layer_j) in combinations(layers, 2):
        print("Comparing layer: {layer_i:s} with {layer_j:s}".format(layer_i=layer_i, layer_j=layer_j))

        Gi = data_cell['graphs'][layer_i]
        Gj = data_cell['graphs'][layer_j]
        dfPCA_i = data_cell['modules-svd']['dfs'][layer_i]
        dfPCA_j = data_cell['modules-svd']['dfs'][layer_j]
        modules_i = data_cell['modules-svd']['modules'][layer_i]
        modules_j = data_cell['modules-svd']['modules'][layer_j]

        for module_i, module_j in [(a, b) for a in modules_i for b in modules_j]:

            id_i = str(module_i['id'])
            id_j = str(module_j['id'])
            name_i = module_i['name']
            name_j = module_j['name']
            print("Comparing: {layer_i:s}-{id_i:s}-{name_i:s} with {layer_j:s}-{id_j:s}-{name_j:s}".format(
                id_i=id_i, id_j=id_j, layer_i=layer_i, layer_j=layer_j, name_i=name_i, name_j=name_j)
            )

            icx = "{:d}c".format(module_i['xy-coords']['x-comp'])
            icy = "{:d}c".format(module_i['xy-coords']['y-comp'])
            icxl, icxh = module_i['xy-coords']['x-values']
            icyl, icyh = module_i['xy-coords']['y-values']

            dfPCA_i_tmp = dfPCA_i.loc[
                (
                    (dfPCA_i[icx] >= icxl) & (dfPCA_i[icx] <= icxh) &
                    (dfPCA_i[icy] >= icyl) & (dfPCA_i[icy] <= icyh)
                ), ['gene', icx, icy]]

            jcx = "{:d}c".format(module_j['xy-coords']['x-comp'])
            jcy = "{:d}c".format(module_j['xy-coords']['y-comp'])
            jcxl, jcxh = module_j['xy-coords']['x-values']
            jcyl, jcyh = module_j['xy-coords']['y-values']
            dfPCA_j_tmp = dfPCA_j.loc[
                (
                    (dfPCA_j[jcx] >= jcxl) & (dfPCA_j[jcx] <= jcxh) &
                    (dfPCA_j[jcy] >= jcyl) & (dfPCA_j[jcy] <= jcyh)
                ), ['gene', jcx, jcy]]

            genes_i = dfPCA_i_tmp.index.to_list()
            genes_j = dfPCA_j_tmp.index.to_list()

            genes_ij = genes_i + genes_j

            # Only genes in this modules
            Gtmp = nx.subgraph(G, genes_ij).copy()

            # Remove intra edges
            remove_intra_edges = [(i, j) for i, j, d in Gtmp.edges(data=True) if d.get('type', None) == 'intra']
            Gtmp.remove_edges_from(remove_intra_edges)

            # Remove isolates
            remove_isolates_nodes = list(nx.isolates(Gtmp))
            Gtmp.remove_nodes_from(remove_isolates_nodes)

            # Jaccard Proximity
            a = set(genes_i)
            b = set(genes_j)
            a_union_b = a.union(b)
            a_inter_b = set(Gtmp.nodes())

            dist = len(a_inter_b) / len(a_union_b)
            r.append((layer_i, layer_j, id_i, id_j, name_i, name_j, dist))

    dfR = pd.DataFrame(r, columns=['layer-i', 'layer-j', 'id-i', 'id-j', 'name-i', 'name-j', 'proximity'])

    ##
    # Export
    ##
    print('Exporting')
    wCSVfile = 'results/proximity/svd-proximity-{celltype:s}-{network:s}-{threshold:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str)
    ensurePathExists(wCSVfile)
    dfR.to_csv(wCSVfile)

if __name__ == '__main__':

    celltype = 'spermatocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5

    # calc_module_proximity(celltype=celltype, network=network, threshold=threshold)

    for celltype in ['spermatocyte', 'enterocyte']:
        calc_module_proximity(celltype=celltype, network=network, threshold=threshold)
