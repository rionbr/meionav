# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM), modules and page rank results and calculates a null model of module page rank.
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

        data_cell = data_cells[celltype]

        for layer in ['HS', 'MM', 'DM']:

            print('Isolate {layer:s} Layer'.format(layer=layer))
            Gt = get_network_layer(G, layer=layer)

            print('Loading SVD for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
            rSVDFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            dfSVD = pd.read_csv(rSVDFile, index_col=0)

            print('Loading PageRank for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
            rPRFile = 'results/pagerank/{celltype:s}/pagerank-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            dfPR = pd.read_csv(rPRFile, index_col=0)

            # Modules
            modules = data_cell['modules-svd']['modules'][layer]

            results = []

            for module in modules:

                id = str(module['id'])
                name = module['name']

                print("Module: {layer:s}-{id:s}-{name:s}".format(id=id, layer=layer, name=name))

                cx = "{:d}c".format(module['xy-coords']['x-comp'])
                cy = "{:d}c".format(module['xy-coords']['y-comp'])
                cxl, cxh = module['xy-coords']['x-values']
                cyl, cyh = module['xy-coords']['y-values']

                # Module df
                dfSVDmod = dfSVD.loc[
                    (
                        (dfSVD[cx] >= cxl) & (dfSVD[cx] <= cxh) &
                        (dfSVD[cy] >= cyl) & (dfSVD[cy] <= cyh)
                    ), ['gene', cx, cy]]
                # Gene List
                genes = dfSVDmod.index

                # Module size
                size = len(dfSVDmod)

                dfPRmod = dfPR.loc[genes, 'page_rank']
                dict_results = dfPRmod.apply(['mean', 'std']).to_dict()
                results.append((celltype, layer, id, name, size, 'real', dfPRmod.tolist()))
                
                for run in range(1, 101):
                    dfPRmodnull = dfPR.sample(n=size)['page_rank']
                    results.append((celltype, layer, id, name, size, run, dfPRmodnull.tolist()))


            dfR = pd.DataFrame(results, columns=['celtype', 'layer', 'mod-id', 'mod-name', 'mod-size', 'run', 'values'])
            ##
            # Export
            ##
            print('Saving results to .CSV')
            wCSVFile = 'results/module_null/{celltype:s}/module-null-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            ensurePathExists(wCSVFile)
            dfR.to_csv(wCSVFile)

    print('Done.')