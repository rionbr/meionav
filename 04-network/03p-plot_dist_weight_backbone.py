# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots networtk weight distribution for all three species
#
# Instructions:
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from utils import get_network_layer


if __name__ == '__main__':

    network = 'complete'  # 'complete'
    threshold = 0.9

    print('Loading Complete Network')
    rGfile_gpickle = 'results/net_{network:s}_mlayer_backbone_modules.gpickle'.format(network='complete')
    G = nx.read_gpickle(rGfile_gpickle)

    #Keep only the backbone
    print('Removing all edges except of backbone edges')
    print('Number of edges: {:,d}'.format(G.number_of_edges()))
    print('Computing edges to remove')
    edges_to_remove = [(i, j) for i, j, d in G.edges(data=True) if (
        ((d.get('type') == 'intra') and (d.get('metric-backbone', False) == False))
    )]
    print('Removing {:,d} edges'.format(len(edges_to_remove)))
    G.remove_edges_from(edges_to_remove)
    print('Number of edges in backbone: {:,d}'.format(G.number_of_edges()))

    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    print('Get edge weights')
    weight_HS = sorted([d['weight'] for i, j, d in HSG.edges(data=True)], reverse=True)
    weight_MM = sorted([d['weight'] for i, j, d in MMG.edges(data=True)], reverse=True)
    weight_DM = sorted([d['weight'] for i, j, d in DMG.edges(data=True)], reverse=True)

    print('Plot')
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # Title
    ax.set_title('Network backbone weight distribution')

    # Plots
    phs, = ax.plot(weight_HS, lw=0, marker='o', ms=4, color='#1f77b4', rasterized=True)
    pmm, = ax.plot(weight_MM, lw=0, marker='o', ms=4, color='#ff7f0e', rasterized=True)
    pdm, = ax.plot(weight_DM, lw=0, marker='o', ms=4, color='#2ca02c', rasterized=True)

    # Insert Distribution Plot
    axin = inset_axes(ax, width='50%', height='50%', loc='lower left', bbox_to_anchor=(.075, .35, .9, .8), bbox_transform=ax.transAxes)
    axin.hist(weight_HS, bins=23, density=True, edgecolor='#1f77b4', facecolor=(0, 0, 0, 0), lw=1)
    axin.hist(weight_MM, bins=23, density=True, edgecolor='#ff7f0e', facecolor=(0, 0, 0, 0), lw=1)
    axin.hist(weight_DM, bins=23, density=True, edgecolor='#2ca02c', facecolor=(0, 0, 0, 0), lw=1)

    # Threshold horizontal line
    ax.axhline(y=threshold, color='red')
    axin.axvline(x=threshold, color='red')

    ax.set_ylabel('Weight')
    ax.set_xlabel('Edge rank')
    ax.set_xscale('log')

    # Legend
    ax.legend(
        handles=(phs, pmm, pdm),
        labels=('Human sapiens', 'Mus musculus', 'Drosophila melanogaster'),
        loc='lower left'
    )

    # Grid
    ax.grid()
    axin.grid()

    plt.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
    file = 'images/img-net-complete-weight-backbone-dist.pdf'
    fig.savefig(file)
