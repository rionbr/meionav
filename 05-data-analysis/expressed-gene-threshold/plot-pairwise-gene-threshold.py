# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots a heatmap of the pairwise jaccard comparison for conserved genes
#
# Instructions:
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
import matplotlib.gridspec as gridspec
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from itertools import combinations
from utils import ensurePathExists


if __name__ == '__main__':

    dfc = pd.read_csv('results/pairwise-celltypes-threshold-similarity.csv.gz', index_col=0)
    dfg = pd.read_csv('results/similarity/similarity-genome.csv', index_col=['layer-i', 'layer-j'])
    # Rename values
    dfc['celltype'] = dfc['celltype'].replace({'enterocyte': 'intestine'})
    # Order by Category
    dfc['celltype'] = pd.Categorical(dfc['celltype'], categories=['spermatogonia', 'spermatocyte', 'spermatid', 'intestine', 'neuron', 'muscle'], ordered=True)

    specie_names = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}

    dfp = dfc.pivot_table(index=['specie_i', 'specie_j', 'celltype'], columns='threshold', values='similarity')

    # Plot
    fig = plt.figure(figsize=(5.4, 8))
    gs = gridspec.GridSpec(ncols=23, nrows=20, figure=fig)
    norm = mpl.colors.Normalize(vmin=0, vmax=1.0)

    cmap = 'seismic'
    #cmap = mpl.colors.LinearSegmentedColormap.from_list(name='br', colors=['#1f77b4', 'white', '#d62728'], N=256)
    ax_y_start = 0
    ax_rx_start, ax_rx_end = 3, 19

    idx = pd.IndexSlice
    #
    #
    #
    ax_y_end = ax_y_start + 1
    axc = fig.add_subplot(gs[ax_y_start: ax_y_end, ax_rx_start: ax_rx_end])
    axc.set_title("Minimum expression threshold", fontsize='medium')
    axc.set_xlabel(r'TPM')
    polygon = plt.Polygon([(0, 0), (10, 0), (10, 1)], color='black')
    axc.add_patch(polygon)
    axc.set_yticks([])
    axc.set_yticklabels([])
    axc.set_xticks(np.arange(0.5, 9.5, 1))
    axc.set_xticklabels(['.01', '.05', '.1', '.2', '.5', '1', '2', '5', '10'])
    axc.spines['top'].set_visible(False)
    axc.spines['bottom'].set_visible(False)
    axc.spines['right'].set_visible(False)
    axc.spines['left'].set_visible(False)
    axc.set_xlim(0, 9)
    ax_y_start = ax_y_end + 1
    #
    #
    #
    print("Germ cell Comparison")
    for specie_i, specie_j in combinations(['HS', 'MM', 'DM'], 2):
        print("Heatmap of {specie_i:s} vs {specie_j:s}".format(specie_i=specie_i, specie_j=specie_j))
        ax_y_block_start = ax_y_block_end = ax_y_start
        #
        # Germ cell Comparison
        #
        tissues = ['spermatogonia', 'spermatocyte', 'spermatid']
        dfpt = dfp.loc[idx[specie_i, specie_j, tissues], :].copy()
        dfpt.sort_index(ascending=True, inplace=True)
        #
        gsim = dfg.loc[(specie_i, specie_j), 'similarity']
        #
        name_i = specie_names[specie_i]
        name_j = specie_names[specie_j]

        ax_y_end = ax_y_start + len(dfpt)
        ax_y_block_end += len(dfpt)
        axc = fig.add_subplot(gs[ax_y_start:ax_y_end, ax_rx_start:ax_rx_end])

        im = axc.imshow(dfpt.values, cmap=cmap, norm=norm, interpolation='nearest', aspect='auto')

        axc.set_xticks(range(len(dfpt.columns)))
        axc.set_xticklabels([])
        axc.set_yticks([0, 1, 2])
        yticklabels = [s.title() for s in dfpt.index.get_level_values(level=2)]
        axc.set_yticklabels(yticklabels)
        axc.yaxis.tick_right()
        #
        ax_y_start += len(dfpt) + 0
        #
        # Other Tissues
        #
        tissues = ['intestine', 'neuron', 'muscle']
        dfpt = dfp.loc[idx[specie_i, specie_j, tissues], :].copy()
        dfpt.sort_index(ascending=True, inplace=True)
        #
        #
        name_i = specie_names[specie_i]
        name_j = specie_names[specie_j]

        ax_y_end = ax_y_start + len(dfpt)
        ax_y_block_end += len(dfpt)
        axo = fig.add_subplot(gs[ax_y_start:ax_y_end, ax_rx_start:ax_rx_end])

        im = axo.imshow(dfpt.values, cmap=cmap, norm=norm, interpolation='nearest', aspect='auto')

        axo.set_xticks(range(len(dfpt.columns)))
        axo.set_xticklabels([])
        axo.set_yticks([0, 1, 2])
        yticklabels = [s.title() for s in dfpt.index.get_level_values(level=2)]
        axo.set_yticklabels(yticklabels)
        axo.yaxis.tick_right()
        #
        #
        #
        axg = fig.add_subplot(gs[ax_y_block_start:ax_y_block_end, 1:3])
        axg.imshow([[gsim]], cmap=cmap, norm=norm, interpolation='nearest', aspect='auto')
        axg.text(x=0.1, y=0, s='{sim:.2f}'.format(sim=gsim), rotation=90, va='center', ha='center')
        axg.set_xticks([0])
        axg.set_yticks([0])
        axg.set_xticklabels([])
        axg.set_yticklabels([])
        if (specie_i == 'HS' and specie_j == 'MM'):
            axg.set_title('Whole genome', fontsize='medium', rotation=90)
        #
        ylabel = "{name_i:s}\nx\n{name_j:s}".format(name_i=name_i, name_j=name_j)
        axg.set_ylabel(ylabel, rotation=0, ha='center', va='center', labelpad=20)
        #
        ax_y_start += len(dfpt) + 0

    #
    # Colorbar    
    #
    cax = fig.add_axes([0.07, 0.04, 0.20, 0.013])
    cbar = fig.colorbar(im, cax=cax, ticks=[0, 0.5, 1.0], orientation='horizontal')
    cax.set_title('Similarity index range', fontsize='small')
    cax.tick_params(labelsize='small')

    #
    # Save
    #
    plt.subplots_adjust(left=0.10, right=0.90, bottom=0.09, top=0.94, wspace=0.7, hspace=0.7)
    wIMGfile = 'images/img-pairwise-conserved-heatmap.pdf'
    ensurePathExists(wIMGfile)
    fig.savefig(wIMGfile)
    plt.close()
