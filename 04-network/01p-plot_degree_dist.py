# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots networtk degree distribution for all three species
#
# Instructions:
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
from utils import get_network_layer, get_network_largest_connected_component
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte

    print('Loading Complete Network')
    rGfile_gpickle = 'results/net-{celltype:s}.gpickle'.format(celltype=celltype)
    G = nx.read_gpickle(rGfile_gpickle)

    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    print('Largest Components')
    HSGc = get_network_largest_connected_component(HSG)
    MMGc = get_network_largest_connected_component(MMG)
    DMGc = get_network_largest_connected_component(DMG)

    print('Get degree')
    degree_HS = sorted([deg for i, deg in HSGc.degree()], reverse=True)
    degree_MM = sorted([deg for i, deg in MMGc.degree()], reverse=True)
    degree_DM = sorted([deg for i, deg in DMGc.degree()], reverse=True)

    print('Plot')
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # Title
    ax.set_title('{celltype:s} net degree dist. (largest conn. comp.)'.format(celltype=celltype.title()))

    # Plots
    phs, = ax.plot(degree_HS, lw=0, marker='o', ms=4, color='#2ca02c', rasterized=True)
    pmm, = ax.plot(degree_MM, lw=0, marker='o', ms=4, color='#7f7f7f', rasterized=True)
    pdm, = ax.plot(degree_DM, lw=0, marker='o', ms=4, color='#ff7f0e', rasterized=True)

    ax.set_ylabel('Number of neighbors')
    ax.set_xlabel('Node rank')
    ax.set_xscale('log')

    # Legend
    ax.legend(
        handles=(phs, pmm, pdm),
        labels=('Human', 'Mouse', 'Insect'),
        loc='upper right'
    )

    # Grid
    ax.grid()

    plt.subplots_adjust(left=0.14, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
    file = 'images/img-net-{celltype:s}-degree-dist.pdf'.format(celltype=celltype)
    fig.savefig(file)
