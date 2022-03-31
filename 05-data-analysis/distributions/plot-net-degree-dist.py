# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots networtk degree distribution for all three species
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
from collections import Counter
from scipy import stats
from utils import get_network_layer, get_network_largest_connected_component, ensurePathExists
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte.")
    parser.add_argument("--network", default='thr', type=str, choices=['full', 'thr'], help="Network to load. Defaults to 'full'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')

    print('Loading {network:s} Network'.format(network=network))
    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
    if network == 'full':
        rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}.gpickle'.format(celltype=celltype, network=network)
    elif network == 'thr':
        rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
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
    layer_data = {
        'HS': {
            'name': 'Human',
            'degree': degree_HS,
            'facecolor': '#2ca02c',
            'edgecolor': '#98df8a',
        },
        'MM': {
            'name': 'Mouse',
            'degree': degree_MM,
            'facecolor': '#7f7f7f',
            'edgecolor': '#c7c7c7',
        },
        'DM': {
            'name': 'Fruit fly',
            'degree': degree_DM,
            'facecolor': '#ff7f0e',
            'edgecolor': '#ffbb78'
        }
    }

    for layer in ['HS', 'MM', 'DM']:

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.1, 4.4))

        #
        name = layer_data[layer]['name']
        degree = layer_data[layer]['degree']
        facecolor = layer_data[layer]['facecolor']
        edgecolor = layer_data[layer]['edgecolor']
        #
        print(layer)
        deg_mean = np.mean(degree)
        deg_std = np.std(degree)
        #
        counter_deg = Counter(degree)
        deg, cnt = zip(*counter_deg.items())
        deg = np.array(deg)
        cnt = np.array(cnt)
        prob = np.array(cnt) / sum(cnt)
        #
        #weights = np.ones(shape=len(degree))
        #hist, bin_edges = np.histogram(degree, bins=24, weights=weights)
        #
        max_degree = max(degree)
        #
        rv_x = np.arange(0, max_degree, 1)
        rv_pmf = stats.norm.pdf(rv_x, loc=deg_mean, scale=deg_std)

        phs, = ax.plot(deg, prob, lw=0, marker='o', ms=4, label=name, color=facecolor, rasterized=False)
        #phs = ax.hist(degree, bins=25, weights=weights_HS)
        rvphs, = ax.plot(rv_x, rv_pmf, lw=2, color='#e377c2', rasterized=False)

        ax.set_title('Degree distribution')
        ax.set_ylabel('P(degree)')
        ax.set_xlabel('Degree rank')
        ax.set_xscale('log')
        #ax.set_yscale('log')

        #ax.grid()
        ax.legend(loc='upper right', handletextpad=0)

        img_path = 'images/degree-dist/{celltype:s}/'.format(celltype=celltype)
        file = img_path + 'img-net-{celltype:s}-{network:s}-degree-dist-{layer:s}.pdf'.format(celltype=celltype, network=network, layer=layer)
        #plt.tight_layout()
        plt.subplots_adjust(left=0.24, right=0.96, bottom=0.16, top=0.86, wspace=0.0)
        ensurePathExists(file)
        fig.savefig(file)
        plt.close()
