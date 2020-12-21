# coding=utf-8
# Author: Rion B Correia
# Date: Aug 20, 2020
#
# Description: Plots Venn Diagrams of the relation of conserved genes
#
# Instructions:
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
#import matplotlib as mpl
#mpl.rcParams['font.family'] = 'Helvetica'
#mpl.rcParams['mathtext.fontset'] = 'cm'
#mpl.rcParams['mathtext.rm'] = 'serif'
#import matplotlib.pyplot as plt
#from matplotlib_venn import venn2
from utils import get_network_layer, ensurePathExists
from itertools import combinations


def venn_count(named_sets):
    """ From: https://stackoverflow.com/questions/15553728/counting-intersections-for-all-combinations-in-a-list-of-sets """
    names = set(named_sets)
    for i in range(1, len(named_sets) + 1):
        for to_intersect in combinations(sorted(named_sets), i):
            others = names.difference(to_intersect)
            intersected = set.intersection(*(named_sets[k] for k in to_intersect))
            unioned = set.union(*(named_sets[k] for k in others)) if others else set()
            yield to_intersect, others, len(intersected - unioned)


if __name__ == '__main__':

    celltypes = ['spermatocyte', 'enterocyte', 'neuron', 'muscle']
    layers = ['HS', 'MM', 'DM']
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    #
    # Conserved
    #
    network = 'conserved'
    data = {
        'HS': {},
        'MM': {},
        'DM': {}
    }
    print('-- Conserved --')
    for celltype in celltypes:
        print('Loading {celltype:s} {network:s} {threshold:s}'.format(celltype=celltype, network='conserved', threshold=threshold_str))

        path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
        rGc_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
        Gc = nx.read_gpickle(rGc_file_gpickle)

        for layer in layers:
            print('Separating layer {layer:s}'.format(layer=layer))
            Gcl = get_network_layer(Gc, layer)
            conserved_genes = set(Gcl.nodes())

            data[layer][celltype] = conserved_genes

    for layer in layers:
        print('Calculating venn {layer:s} conserved'.format(layer=layer))

        ns = data[layer]

        for intersected, unioned, count in venn_count(ns):
            print('|{}{}| = {}'.format(' & '.join(sorted(intersected)), ' - ' + ' - '.join(sorted(unioned)) if unioned else '', count))

    #
    # Non-conserved
    #
    print('-- Non-Conserved --')
    for celltype in celltypes:
        print('Loading {celltype:s} (non-){network:s} {threshold:s}'.format(celltype=celltype, network='conserved', threshold=threshold_str))

        path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
        rGt_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='thr', threshold=threshold_str)
        rGc_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
        Gt = nx.read_gpickle(rGt_file_gpickle)
        Gc = nx.read_gpickle(rGc_file_gpickle)

        for layer in layers:
            print('Separating layer {layer:s}'.format(layer=layer))
            Gtl = get_network_layer(Gt, layer)
            Gcl = get_network_layer(Gc, layer)
            
            conserved_genes = set(Gcl.nodes())
            threshold_genes = set(Gtl.nodes())
            non_conserved_genes = threshold_genes.difference(conserved_genes)

            data[layer][celltype] = non_conserved_genes

    for layer in layers:
        print('Calculating venn {layer:s} non-conserved'.format(layer=layer))

        ns = data[layer]

        for intersected, unioned, count in venn_count(ns):
            print('|{}{}| = {}'.format(' & '.join(sorted(intersected)), ' - ' + ' - '.join(sorted(unioned)) if unioned else '', count))
