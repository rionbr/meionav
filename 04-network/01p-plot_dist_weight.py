# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots networtk weight distribution for all three species
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
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    threshold = args.threshold

    print('Loading Complete Network')
    rGfile_gpickle = 'results/net-{celltype:s}.gpickle'.format(celltype=celltype)
    G = nx.read_gpickle(rGfile_gpickle)

    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    print('Get edge weights')
    weight_HS = sorted([d['weight'] for i, j, d in HSG.edges(data=True)], reverse=True)
    weight_MM = sorted([d['weight'] for i, j, d in MMG.edges(data=True)], reverse=True)
    weight_DM = sorted([d['weight'] for i, j, d in DMG.edges(data=True)], reverse=True)

    hist_weights_HS = np.ones_like(weight_HS) / len(weight_HS)
    hist_weights_MM = np.ones_like(weight_MM) / len(weight_MM)
    hist_weights_DM = np.ones_like(weight_DM) / len(weight_DM)

    print('Plot')
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # Title
    ax.set_title('{celltype:s} network weight dist.'.format(celltype=celltype.title()))

    # Plots
    phs, = ax.plot(weight_HS, lw=0, marker='o', ms=4, color='#2ca02c', rasterized=True)
    pmm, = ax.plot(weight_MM, lw=0, marker='o', ms=4, color='#7f7f7f', rasterized=True)
    pdm, = ax.plot(weight_DM, lw=0, marker='o', ms=4, color='#ff7f0e', rasterized=True)

    # Insert Distribution Plot
    axin = inset_axes(ax, width='50%', height='50%', loc='lower left', bbox_to_anchor=(.11, .325, .9, .8), bbox_transform=ax.transAxes)
    axin.hist(weight_HS, bins=20, weights=hist_weights_HS, edgecolor='#2ca02c', facecolor=(0, 0, 0, 0), lw=1)
    axin.hist(weight_MM, bins=20, weights=hist_weights_MM, edgecolor='#7f7f7f', facecolor=(0, 0, 0, 0), lw=1)
    axin.hist(weight_DM, bins=20, weights=hist_weights_DM, edgecolor='#ff7f0e', facecolor=(0, 0, 0, 0), lw=1)

    # Threshold horizontal line
    ax.axhline(y=threshold, color='red')
    axin.axvline(x=threshold, color='red')

    ax.set_ylabel('Weight')
    ax.set_xlabel('Edge rank')
    ax.set_xscale('log')

    axin.set_ylabel('Probability', fontsize='small')
    axin.set_xlabel('Weight', fontsize='small')
    axin.set_xticks([0.2, 0.5, 1.0])

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
    file = 'images/img-net-{celltype:s}-weight-dist.pdf'.format(celltype=celltype)
    fig.savefig(file)
    plt.close()