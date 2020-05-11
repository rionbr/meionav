# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads the SVD results and calculats a similarity between SVD Modules across layers.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
from data import *


if __name__ == '__main__':

    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    layers = ['HS', 'MM', 'DM']
    celltype_i = 'enterocyte'
    celltype_j = 'spermatocyte'
    data_i = data_cells[celltype_i]
    data_j = data_cells[celltype_j]

    for layer in layers:

        print('Loading SVD for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype_i, network=network, threshold=threshold_str, layer=layer))
        rPCAeFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype_i, network=network, threshold=threshold_str, layer=layer)
        dfPCAe = pd.read_csv(rPCAeFile, index_col=0, encoding='utf-8')
        data_i['modules-svd']['files'][layer] = rPCAeFile
        data_i['modules-svd']['dfs'][layer] = dfPCAe

        print('Loading SVD for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype_j, network=network, threshold=threshold_str, layer=layer))
        rPCAsFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype_j, network=network, threshold=threshold_str, layer=layer)
        dfPCAs = pd.read_csv(rPCAsFile, index_col=0, encoding='utf-8')
        data_j['modules-svd']['files'][layer] = rPCAsFile
        data_j['modules-svd']['dfs'][layer] = dfPCAs

    r = []
    for layer in layers:
        print("Comparing {celltype_i:s} with {celltype_j:s} of layer: {layer:s}".format(celltype_i=celltype_i, celltype_j=celltype_j, layer=layer))

        Gi = data_i['graphs'][layer]
        Gj = data_j['graphs'][layer]
        dfPCA_i = data_i['modules-svd']['dfs'][layer]
        dfPCA_j = data_j['modules-svd']['dfs'][layer]
        modules_i = data_i['modules-svd']['modules'][layer]
        modules_j = data_j['modules-svd']['modules'][layer]

        for module_i, module_j in [(a, b) for a in modules_i for b in modules_j]:

            id_i = str(module_i['id'])
            id_j = str(module_j['id'])
            name_i = module_i['name']
            name_j = module_j['name']
            print("Comparing: {layer:s}-{id_i:s}-{name_i:s} with {layer:s}-{id_j:s}-{name_j:s}".format(
                id_i=id_i, id_j=id_j, layer=layer, name_i=name_i, name_j=name_j)
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

            # Jaccard Proximity
            a = set(genes_i)
            b = set(genes_j)
            a_union_b = a.union(b)
            a_inter_b = a.intersection(b)
            dist = len(a_inter_b) / len(a_union_b)
            r.append((layer, celltype_i, celltype_j, id_i, id_j, name_i, name_j, dist))

    dfR = pd.DataFrame(r, columns=['layer', 'celltype-i', 'celltype-j', 'id-i', 'id-j', 'name-i', 'name-j', 'proximity'])

    ##
    # Export
    ##
    print('Exporting')
    wCSVfile = 'results/proximity/svd-proximity-layers-{network:s}-{threshold:s}.csv.gz'.format(network=network, threshold=threshold_str)
    ensurePathExists(wCSVfile)
    dfR.to_csv(wCSVfile)
