# coding=utf-8
# Author: Rion B Correia
# Date: March 30, 2021
#
# Description: Calculates networks backbone stats
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from distanceclosure.utils import dist2prox


if __name__ == '__main__':

    celltypes = ['spermatocyte', 'enterocyte']
    layers = ['HS', 'MM', 'DM']
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    r = []
    for celltype in celltypes:

        print("Celltype: {celltype:s}".format(celltype=celltype))
        path = "../../04-network/results/network-closure/{celltype:s}/".format(celltype=celltype)

        for layer in layers:

            print("Layer: {layer:s}".format(layer=layer))
            rGfile = path + "net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle".format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            G = nx.read_gpickle(rGfile)

            nodes = G.number_of_nodes()
            edges = G.number_of_edges()
            density = nx.density(G)
            are_metric = len([(i, j) for i, j, d in G.edges(data=True) if d.get('is_metric') is True])
            are_ultrametric = len([(i, j) for i, j, d in G.edges(data=True) if d.get('is_ultrametric') is True])

            edge_distortion = [abs(dist2prox(d['metric_distance']) - dist2prox(d['distance'])) for i, j, d in G.edges(data=True)]
            total_distortion = sum(edge_distortion)
            norm_distortion = (2 * total_distortion) / (nodes * (nodes - 1))

            r.append([
                celltype,
                network,
                threshold,
                layer,
                nodes,
                edges,
                density,
                are_metric,
                are_ultrametric,
                total_distortion,
                norm_distortion,
            ])

    cols = ['celltype', 'network', 'threshold', 'layer', 'n-nodes', 'n-edges', 'density', 'n-edges-metric', 'n-edges-ultrametric', 'distortion', 'distortion-norm']
    df = pd.DataFrame(r, columns=cols)

    #Export
    df.to_csv("results/net-backbone-stats.csv")