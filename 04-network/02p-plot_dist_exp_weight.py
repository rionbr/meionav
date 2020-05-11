# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots network (experimental evidence) weight distribution for all three species
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
from utils import get_network_layer
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')


    #
    print('Loading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
    rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    print('Get edge weights')
    exp_weight_HS = sorted([d['experiments'] / 1000 for i, j, d in HSG.edges(data=True) if 'experiments' in d], reverse=True)
    exp_weight_MM = sorted([d['experiments'] / 1000 for i, j, d in MMG.edges(data=True) if 'experiments' in d], reverse=True)
    exp_weight_DM = sorted([d['experiments'] / 1000 for i, j, d in DMG.edges(data=True) if 'experiments' in d], reverse=True)

    print('Plot')
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # Title
    ax.set_title("{celltype:s} {network:s}-{threshold:.1f} net. 'experiment' weight dist.".format(celltype=celltype.title(), network=network, threshold=threshold))

    # Plots
    phs, = ax.plot(exp_weight_HS, lw=0, marker='o', ms=4, color='#2ca02c', rasterized=True)
    pmm, = ax.plot(exp_weight_MM, lw=0, marker='o', ms=4, color='#7f7f7f', rasterized=True)
    pdm, = ax.plot(exp_weight_DM, lw=0, marker='o', ms=4, color='#ff7f0e', rasterized=True)

    # Insert Distribution Plot
    hist_weights_HS = np.ones_like(exp_weight_HS) / len(exp_weight_HS)
    hist_weights_MM = np.ones_like(exp_weight_MM) / len(exp_weight_MM)
    hist_weights_DM = np.ones_like(exp_weight_DM) / len(exp_weight_DM)
    axin = inset_axes(ax, width='42%', height='42%', loc='lower left', bbox_to_anchor=(.11, .325, .9, .8), bbox_transform=ax.transAxes)
    axin.hist(exp_weight_HS, bins=23, weights=hist_weights_HS, edgecolor='#2ca02c', facecolor=(0, 0, 0, 0), lw=1)
    axin.hist(exp_weight_MM, bins=23, weights=hist_weights_MM, edgecolor='#7f7f7f', facecolor=(0, 0, 0, 0), lw=1)
    axin.hist(exp_weight_DM, bins=23, weights=hist_weights_DM, edgecolor='#ff7f0e', facecolor=(0, 0, 0, 0), lw=1)

    ax.set_ylabel('Weight')
    ax.set_xlabel('Edge rank')
    ax.set_xscale('log')

    axin.set_ylabel('Probability', fontsize='small')
    axin.set_xlabel('Weight', fontsize='small')

    # Legend
    ax.legend(
        handles=(phs, pmm, pdm),
        labels=('Human', 'Mouse', 'Insect'),
        loc='lower left'
    )

    # Grid
    ax.grid()
    axin.grid()

    plt.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
    file = 'images/img-net-{celltype:s}-{network:s}-{threshold:s}-exp-weight-dist.pdf'.format(celltype=celltype, network=network, threshold=threshold_str)
    fig.savefig(file)
