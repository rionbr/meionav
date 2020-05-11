# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots a graph that shows different network weight thresholding and the loss of core genes
#
# Instructions:
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
from utils import get_network_layer
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte

    print('Reading Network')
    rGfile_gpickle = 'results/net-{celltype:s}.gpickle'.format(celltype=celltype)
    G = nx.read_gpickle(rGfile_gpickle)

    DMG = get_network_layer(G, 'DM')

    core_DM = nx.get_node_attributes(DMG, name='core')
    gene_DM = nx.get_node_attributes(DMG, name='label')
    df_DM_m = pd.DataFrame(data={'gene': gene_DM, 'core': core_DM})
    df_DM_m['core'] = df_DM_m['core'].fillna(False)

    r = []
    Gt = DMG.copy()
    r.append([
        None,
        Gt.number_of_nodes(),
        Gt.number_of_edges(),
        len([i for i, d in Gt.nodes(data=True) if d.get('core', False) == True])
    ])
    for threshold in np.linspace(0.0, 1, 20, False):
        print('Threshold: {:.2f}'.format(threshold))

        # Remove Edges
        edges_to_remove = [(i, j) for i, j, d in Gt.edges(data=True) if d['weight'] < threshold]
        Gt.remove_edges_from(edges_to_remove)

        # Remove Isolates
        nodes_to_remove = list(nx.isolates(Gt))
        Gt.remove_nodes_from(nodes_to_remove)
        #
        r.append([
            threshold,
            Gt.number_of_nodes(),
            Gt.number_of_edges(),
            len([i for i, d in Gt.nodes(data=True) if d.get('core', False) == True])
        ])

    dfR = pd.DataFrame(r, columns=['threshold', 'n-nodes', 'n-edges', 'n-core'])

    print("Plot")
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    pdm, = ax.plot(dfR['threshold'], dfR['n-core'], color='#ff7f0e', marker='o', lw=2)
    ax.set_title('{celltype:s} network threshold. Core gene loss'.format(celltype=celltype.title()))
    ax.set_ylabel('Number of core genes')
    ax.set_xlabel('Edge weight threshold')
    ax.grid()
    ax.set_xticks(np.linspace(0, 1, 5))
    #ax.set_xticklabels(np.linspace(0, 1, 6))


    # Legend
    ax.legend(
        handles=(pdm,),
        labels=('Insect',),
        loc='lower left'
    )

    plt.subplots_adjust(left=0.12, right=0.96, bottom=0.12, top=0.92, wspace=0, hspace=0)
    file = 'images/img-net-{celltype:s}-thr-cores.pdf'.format(celltype=celltype)
    fig.savefig(file)
    plt.close()
