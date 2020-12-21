# coding=utf-8
# Author: Rion B Correia
# Date: Aug 20, 2020
#
# Description: Plots Venn Diagrams of the relation between threshold and conserved genes
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
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from utils import get_network_layer, ensurePathExists


if __name__ == '__main__':

    celltype = 'enterocyte'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    network = 'thr'
    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
    rG_thr_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    print('Load {celltype:s} {network:s} graph'.format(celltype=celltype, network=network))
    Gt = nx.read_gpickle(rG_thr_file_gpickle)

    print('Separate Layers')
    HSGt = get_network_layer(Gt, 'HS')
    MMGt = get_network_layer(Gt, 'MM')
    DMGt = get_network_layer(Gt, 'DM')

    network = 'conserved'
    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
    rG_con_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    print('Load {celltype:s} {network:s} graph'.format(celltype=celltype, network=network))
    Gc = nx.read_gpickle(rG_con_file_gpickle)

    print('Separate Layers')
    HSGc = get_network_layer(Gc, 'HS')
    MMGc = get_network_layer(Gc, 'MM')
    DMGc = get_network_layer(Gc, 'DM')

    dict_data = {
        'HS': {
            'graph_thr': HSGt,
            'graph_conserved': HSGc,
            'specie': 'Human',
            'facecolor': '#98df8a',
            'edgecolor': '#2ca02c'
        },
        'MM': {
            'graph_thr': MMGt,
            'graph_conserved': MMGc,
            'specie': 'Mouse',
            'facecolor': '#c7c7c7',
            'edgecolor': '#7f7f7f',
        },
        'DM': {
            'graph_thr': DMGt,
            'graph_conserved': DMGc,
            'specie': 'Insect',
            'facecolor': '#ffbb78',
            'edgecolor': '#ff7f0e'
        }
    }

    print('Plot')
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(7, 3.))

    for ax, layer in zip(axes, ['HS', 'MM', 'DM']):

        Gtt = dict_data[layer]['graph_thr']
        Gtc = dict_data[layer]['graph_conserved']
        specie = dict_data[layer]['specie']
        facecolor = dict_data[layer]['facecolor']
        edgecolor = dict_data[layer]['edgecolor']

        set_thr_genes = set(Gtt.nodes())
        set_conserved_genes = set(Gtc.nodes())

        title = "{specie:s} {celltype:s}\n(thr={threshold:.1f})".format(celltype=celltype, specie=specie, threshold=threshold)
        ax.set_title(title)

        def label_formatter(t):
            return '{t:,d}'.format(t=t)
        # Plots
        v = venn2(ax=ax,
                  subsets=[set_thr_genes, set_conserved_genes],
                  set_labels=('Threshold', 'Conserved'),
                  #set_colors=[],
                  alpha=1.0,
                  subset_label_formatter=label_formatter,
                  normalize_to=1.0)
        # Patches
        patch_A = v.get_patch_by_id('10')
        patch_B = v.get_patch_by_id('01')
        patch_AB = v.get_patch_by_id('11')
        #
        patch_A.set_facecolor(facecolor)
        patch_A.set_edgecolor(edgecolor)
        patch_A.set_linestyle('solid')
        #
        patch_B.set_facecolor(facecolor)
        patch_B.set_edgecolor(edgecolor)
        patch_B.set_linestyle('solid')
        #
        patch_AB.set_facecolor('#ff9896')  # red
        patch_AB.set_edgecolor('#d62728')
        patch_AB.set_linestyle('dashed')

        # Label
        label_A = v.get_label_by_id('10')
        label_B = v.get_label_by_id('01')
        label_AB = v.get_label_by_id('11')
        #
        label_B.set_visible(False)
        # Sub Labels
        set_label_A = v.set_labels[0]
        set_label_B = v.set_labels[1]
        #
        #set_label_all.set_horizontalalignment('left')
        #set_label_conserved.set_horizontalalignment('right')
        set_label_A.set_fontsize('medium')
        set_label_B.set_fontsize('medium')
        #

    plt.tight_layout()
    img_path = 'images/venn-conserved/'
    file = img_path + 'img-venn-{celltype:s}-thr-vs-conserved-{threshold:s}.pdf'.format(celltype=celltype, threshold=threshold_str)
    ensurePathExists(file)
    fig.savefig(file)
    plt.close()
