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
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from utils import get_network_layer, ensurePathExists


if __name__ == '__main__':

    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype='spermatocyte')
    rG_spermatocyte_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype='spermatocyte', network='conserved', threshold=threshold_str)

    print('Load Spermatocyte graph')
    Gs = nx.read_gpickle(rG_spermatocyte_file_gpickle)

    print('Separate Layers')
    HSGs = get_network_layer(Gs, 'HS')
    MMGs = get_network_layer(Gs, 'MM')
    DMGs = get_network_layer(Gs, 'DM')

    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype='enterocyte')
    rG_enterocyte_file_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype='enterocyte', network='conserved', threshold=threshold_str)
    print('Load Enterocyte graph')
    Ge = nx.read_gpickle(rG_enterocyte_file_gpickle)

    print('Separate Layers')
    HSGe = get_network_layer(Ge, 'HS')
    MMGe = get_network_layer(Ge, 'MM')
    DMGe = get_network_layer(Ge, 'DM')

    dict_data = {
        'HS': {
            'graph_spermatocyte': HSGs,
            'graph_enterocyte': HSGe,
            'specie': 'Human',
            'facecolor': '#ff9896',
            'edgecolor': '#2ca02c'
        },
        'MM': {
            'graph_spermatocyte': MMGs,
            'graph_enterocyte': MMGe,
            'specie': 'Mouse',
            'facecolor': '#ff9896',
            'edgecolor': '#7f7f7f',
        },
        'DM': {
            'graph_spermatocyte': DMGs,
            'graph_enterocyte': DMGe,
            'specie': 'Insect',
            'facecolor': '#ff9896',
            'edgecolor': '#ff7f0e'
        }
    }

    print('Plot')
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(7, 3.))

    for ax, layer in zip(axes, ['HS', 'MM', 'DM']):

        Gts = dict_data[layer]['graph_spermatocyte']
        Gte = dict_data[layer]['graph_enterocyte']
        specie = dict_data[layer]['specie']
        facecolor = dict_data[layer]['facecolor']
        edgecolor = dict_data[layer]['edgecolor']

        set_spermatocyte_genes = set(Gts.nodes())
        set_enterocyte_genes = set(Gte.nodes())

        title = "{specie:s}\n(conserved, thr={threshold:.1f})".format(specie=specie, threshold=threshold)
        ax.set_title(title)

        def label_formatter(t):
            return '{t:,d}'.format(t=t)
        # Plots
        v = venn2(ax=ax,
                  subsets=[set_spermatocyte_genes, set_enterocyte_genes],
                  set_labels=('Spermatocyte', 'Enterocyte'),
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
        patch_AB.set_facecolor('#ffb6b5')  # red
        patch_AB.set_edgecolor(edgecolor)
        patch_AB.set_linestyle('dashed')

        # Label
        label_A = v.get_label_by_id('10')
        label_B = v.get_label_by_id('01')
        label_AB = v.get_label_by_id('11')
        #
        #label_conserved.set_visible(False)

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
    file = img_path + 'img-venn-conserved-{threshold:s}-spermatocyte-vs-enterocyte.pdf'.format(threshold=threshold_str)
    ensurePathExists(file)
    fig.savefig(file)
    plt.close()
