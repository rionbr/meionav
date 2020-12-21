# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots networtk weight distribution for all three species
#
# Instructions:
#
import math
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
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
import argparse
from utils import get_network_layer, ensurePathExists


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    network = 'full'
    attribute = 'combined_score'
    threshold = args.threshold

    print('Loading Full Network')
    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
    rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}.gpickle'.format(celltype=celltype, network=network)
    G = nx.read_gpickle(rGfile_gpickle)

    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    print('Get edge weights')
    values_HS = sorted([d[attribute] / 1000 for i, j, d in HSG.edges(data=True) if attribute in d], reverse=True)
    values_MM = sorted([d[attribute] / 1000 for i, j, d in MMG.edges(data=True) if attribute in d], reverse=True)
    values_DM = sorted([d[attribute] / 1000 for i, j, d in DMG.edges(data=True) if attribute in d], reverse=True)

    hist_values_HS = np.ones_like(values_HS) / len(values_HS)
    hist_values_MM = np.ones_like(values_MM) / len(values_MM)
    hist_values_DM = np.ones_like(values_DM) / len(values_DM)

    print('Plot')
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # Title
    ax.set_title("{celltype:s} full net 'combined-score' dist.".format(celltype=celltype.title()))

    # Plots
    phs, = ax.plot(np.arange(1, len(values_HS) + 1), values_HS, lw=0, marker='o', ms=4, color='#2ca02c', rasterized=True, zorder=5)
    pmm, = ax.plot(np.arange(1, len(values_MM) + 1), values_MM, lw=0, marker='o', ms=4, color='#7f7f7f', rasterized=True, zorder=4)
    pdm, = ax.plot(np.arange(1, len(values_DM) + 1), values_DM, lw=0, marker='o', ms=4, color='#ff7f0e', rasterized=True, zorder=3)

    max_value = math.ceil(max(values_HS + values_MM + values_DM))
    bins = np.linspace(0, max_value, 23, endpoint=True)

    # Insert Distribution Plot
    axin = inset_axes(ax, width='50%', height='50%', loc='lower left', bbox_to_anchor=(.11, .325, .9, .8), bbox_transform=ax.transAxes)
    axin.hist(values_HS, bins=bins, weights=hist_values_HS, edgecolor='#2ca02c', facecolor=(0, 0, 0, 0), lw=1, zorder=5)
    axin.hist(values_MM, bins=bins, weights=hist_values_MM, edgecolor='#7f7f7f', facecolor=(0, 0, 0, 0), lw=1, zorder=4)
    axin.hist(values_DM, bins=bins, weights=hist_values_DM, edgecolor='#ff7f0e', facecolor=(0, 0, 0, 0), lw=1, zorder=3)

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
    ax.grid(zorder=1)
    axin.grid(zorder=1)

    plt.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
    img_path = 'images/net-edge-attributes/{celltype:s}/'.format(celltype=celltype)
    file = img_path + 'img-net-{celltype:s}-full-edge-{attribute:s}-dist.pdf'.format(celltype=celltype, attribute=attribute)
    ensurePathExists(file)
    fig.savefig(file)
    plt.close()
