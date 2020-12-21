# coding=utf-8
# Author: Rion B Correia
# Date: July 22, 2020
#
# Description: Loads data from 'calc--assortativity' and plots
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def annotate_heatmap(dft):
    for x in [0, 1]:
        for y in [0, 1]:
            plt.text(x, y, '{:.0%}'.format(dft.iloc[x, y]), horizontalalignment='center', verticalalignment='center', fontsize='small')


if __name__ == '__main__':

    celltype = 'spermatocyte'
    layers = ['HS', 'MM', 'DM']

    df = pd.read_csv('results/csv-null-model-assortativity.csv', index_col=0)

    data = {
        'HS': {
            'name': 'Human',
            'facecolor': '#2ca02c',
            'edgecolor': '#98df8a',
        },
        'MM': {
            'name': 'Mouse',
            'facecolor': '#7f7f7f',
            'edgecolor': '#c7c7c7',
        },
        'DM': {
            'name': 'Insect',
            'facecolor': '#ff7f0e',
            'edgecolor': '#ffbb78',
        }
    }

    # Transform into matrix
    for layer in layers:
        dft = df.loc[(df['layer'] == layer), ['con-con', 'con-non', 'non-con', 'non-non']]
        dftm = pd.DataFrame(dft.values.reshape(2, 2), columns=['conserved', 'non-conserved'], index=['conserved', 'non-conserved'])
        dftm = dftm / dftm.sum().sum()
        data[layer]['df'] = dftm

    facecolors = ['#2ca02c', '#7f7f7f', '#ff7f0e']
    edgecolors = []
    measures = ['degree-centrality', 'page-rank']
    titles = ['Degree centrality', 'Page rank']

    # Plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.8, 3.6))
    #
    inds = [1, 2, 3]
    heights = df['assortativity'].tolist()
    bars = ax.bar(inds, heights, width=0.45, zorder=7)
    for patch, facecolor in zip(bars, facecolors):
        patch.set_facecolor(facecolor)
        patch.set_edgecolor('black')

    ax.axhline(0, color='black')
    #
    ax.set_title('Assortativity')
    ax.set_ylabel('Correlation')
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['Human', 'Mouse', 'Insect'])
    ax.set_xlim(0.4, 3.6)
    ax.set_ylim(-0.025, 0.5)
    #
    ticks, ticklabels = [0, 1], ['C', 'N']
    #
    # Insets
    axinhs = inset_axes(ax, width="30%", height="30%", bbox_to_anchor=(0.02, 0.4, 1.0, 1.0), bbox_transform=ax.transAxes, loc=3)

    cmaphs = LinearSegmentedColormap.from_list(colors=['white', '#2ca02c'], name='hs')
    axinhs.imshow(data['HS']['df'], cmap=cmaphs)
    axinhs.set_yticks(ticks)
    axinhs.set_xticks(ticks)
    axinhs.set_xticklabels(ticklabels)
    axinhs.set_yticklabels(ticklabels)
    annotate_heatmap(data['HS']['df'])

    axinmm = inset_axes(ax, width="30%", height="30%", bbox_to_anchor=(0.335, 0.37, 1.0, 1.0), bbox_transform=ax.transAxes, loc=3)

    cmapmm = LinearSegmentedColormap.from_list(colors=['white', '#7f7f7f'], name='mm')
    axinmm.imshow(data['MM']['df'], cmap=cmapmm)
    axinmm.set_yticks(ticks)
    axinmm.set_xticks(ticks)
    axinmm.set_xticklabels(ticklabels)
    axinmm.set_yticklabels(ticklabels)
    annotate_heatmap(data['MM']['df'])

    axindm = inset_axes(ax, width="30%", height="30%", bbox_to_anchor=(0.65, 0.63, 1.0, 1.0), bbox_transform=ax.transAxes, loc=3)

    cmapdm = LinearSegmentedColormap.from_list(colors=['white', '#ff7f0e'], name='dm')
    axindm.imshow(data['DM']['df'], cmap=cmapdm)
    axindm.set_yticks(ticks)
    axindm.set_xticks(ticks)
    axindm.set_xticklabels(ticklabels)
    axindm.set_yticklabels(ticklabels)
    annotate_heatmap(data['DM']['df'])
    #
    #
    ax.grid(axis='y')

    wIMGfile = 'images/network-{network:s}-assort.pdf'.format(network='thr-conserved')
    ensurePathExists(wIMGfile)
    #
    #plt.subplots_adjust(left=0.03, right=0.90, bottom=0.12, top=0.87, wspace=0.0, hspace=0.0)
    #plt.tight_layout()
    plt.savefig(wIMGfile)
    plt.close()

