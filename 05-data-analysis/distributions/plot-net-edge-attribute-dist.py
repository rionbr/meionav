# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots thresholded network edge attribute distribution for all three species
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
#
mpl.rc('font', size=16)  # controls default text sizes
mpl.rc('axes', titlesize=20)  # fontsize of the axes title
mpl.rc('axes', labelsize=16)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=14)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=14)  # fontsize of the tick labels
mpl.rc('legend', fontsize=16)  # legend fontsize
mpl.rc('figure', titlesize=20)  # fontsize of the figure title
#
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

import matplotlib.pyplot as plt
import argparse
from utils import get_network_layer, get_network_by_attribute, ensurePathExists


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    parser.add_argument("--biotype", default=None, type=str, choices=['protein_coding'], help="Filter nodes by biotype (e.g., protein-coding)")
    default_attributes = ['combined_score', 'textmining', 'database', 'experiments', 'coexpression', 'cooccurence', 'fusion']  # , 'neighborhood']
    parser.add_argument("--attribute", default=default_attributes, type=str, help="Which attribute to plot. Defaults to StringDB types")
    #
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    biotype = args.biotype
    attributes = args.attribute if isinstance(args.attribute, list) else [args.attribute]

    print('Loading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
    if network == 'thr':
        rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    if biotype is not None:
        print('Filter by node attribute')
        HSG = get_network_by_attribute(HSG, attribute='biotype', value=biotype)
        MMG = get_network_by_attribute(MMG, attribute='biotype', value=biotype)
        DMG = get_network_by_attribute(DMG, attribute='biotype', value=biotype)

    for attribute in attributes:

        if attribute == 'combined_score':
            attribute_str = 'combined-score'
        else:
            attribute_str = attribute

        title = "{celltype:s} {network:s}-{threshold:.1f} net '{attribute:s}' dist.".format(celltype=celltype.title(), network=network, threshold=threshold, attribute=attribute_str)

        values_HS = sorted([d[attribute] / 1000 for i, j, d in HSG.edges(data=True) if attribute in d], reverse=True)
        values_MM = sorted([d[attribute] / 1000 for i, j, d in MMG.edges(data=True) if attribute in d], reverse=True)
        values_DM = sorted([d[attribute] / 1000 for i, j, d in DMG.edges(data=True) if attribute in d], reverse=True)

        # Plot
        print('Plot {celltype:s}-{network:s}-{threshold:s}-edge-{attribute:s}'.format(celltype=celltype, network=network, threshold=threshold_str, attribute=attribute_str))
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7.5, 6))

        # Title
        ax.set_title(title)

        # Plots
        phs, = ax.plot(np.arange(1, len(values_HS) + 1), values_HS, lw=0, marker='o', ms=6, color='#2ca02c', rasterized=True, zorder=5)
        pmm, = ax.plot(np.arange(1, len(values_MM) + 1), values_MM, lw=0, marker='o', ms=6, color='#7f7f7f', rasterized=True, zorder=4)
        pdm, = ax.plot(np.arange(1, len(values_DM) + 1), values_DM, lw=0, marker='o', ms=6, color='#ff7f0e', rasterized=True, zorder=3)

        max_value = math.ceil(max(values_HS + values_MM + values_DM))
        bins = np.linspace(0, max_value, 23, endpoint=True)

        # Insert Distribution Plot
        hist_values_HS = np.ones_like(values_HS) / len(values_HS)
        hist_values_MM = np.ones_like(values_MM) / len(values_MM)
        hist_values_DM = np.ones_like(values_DM) / len(values_DM)
        axin = inset_axes(ax, width='40%', height='40%', loc='lower left', bbox_to_anchor=(.15, .35, .9, .8), bbox_transform=ax.transAxes)
        axin.hist(values_HS, bins=bins, weights=hist_values_HS, edgecolor='#2ca02c', facecolor=(0, 0, 0, 0), lw=1.2, zorder=5)
        axin.hist(values_MM, bins=bins, weights=hist_values_MM, edgecolor='#7f7f7f', facecolor=(0, 0, 0, 0), lw=1.2, zorder=4)
        axin.hist(values_DM, bins=bins, weights=hist_values_DM, edgecolor='#ff7f0e', facecolor=(0, 0, 0, 0), lw=1.2, zorder=3)

        ax.set_ylabel('Weight')
        ax.set_xlabel('Edge rank')
        ax.set_xscale('log')
        axin.set_ylabel('Probability', fontsize='small')
        axin.set_xlabel('Weight', fontsize='small')

        # Legend
        ax.legend(
            handles=(phs, pmm, pdm),
            labels=('Human', 'Mouse', 'Fruit fly'),
            loc='lower left'
        )

        # Grid
        #ax.grid(zorder=1)
        #axin.grid(zorder=1)

        plt.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
        path_img = 'images/net-edge-attributes/{celltype:s}/'.format(celltype=celltype)
        if network == 'thr':
            file = path_img + 'img-net-{celltype:s}-{network:s}-{threshold:s}-edge-{attribute:s}-dist.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, attribute=attribute_str)
        ensurePathExists(file)
        fig.savefig(file)
        plt.close()
