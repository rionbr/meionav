# coding=utf-8
# Author: Rion B Correia
# Date: Sept 08, 2020
#
# Description: Reads the jaccard results and plots the results.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from utils import ensurePathExists


if __name__ == '__main__':

    #
    print('Loading .csv')
    dfS = pd.read_csv('results/network-jaccard-species.csv', index_col=0)
    dfS['celltype'] = pd.Categorical(dfS['celltype'], categories=['spermatocyte', 'enterocyte', 'neuron', 'muscle'], ordered=True)
    dfS = dfS.set_index(['layer_i', 'layer_j', 'celltype'])
    #
    dfC = pd.read_csv('results/network-jaccard-celltypes.csv', index_col=0)
    dfC['layer'] = pd.Categorical(dfC['layer'], categories=['HS', 'MM', 'DM'], ordered=True)
    dfC = dfC.set_index(['celltype_i', 'celltype_j', 'layer'])



    print('Plot')
    fig = plt.figure(figsize=(3.6, 2.1))
    gs = gridspec.GridSpec(ncols=5, nrows=1, figure=fig, )
    axl = fig.add_subplot(gs[:, 0:3])
    axl.set_anchor('N')
    axr = fig.add_subplot(gs[:, 3:5])
    axr.set_anchor('N')
    cax = fig.add_axes([0.68, 0.11, 0.025, 0.3])

    #
    cmap = cm.get_cmap('seismic')
    #cmap.set_bad(color='#bdbdbd')

    idx = pd.IndexSlice
    # Index Data to plot
    dfl = dfC.loc[idx['spermatocyte', :, :], :].unstack(level=1)
    # reorder columns
    dfl = dfl.reindex(columns=['enterocyte', 'neuron', 'muscle'], level=1)
    print(dfl)
    dfr = dfS.loc[idx['HS', ['MM', 'DM'], 'spermatocyte'], :].T

    iml = axl.imshow(dfl.values, cmap=cmap, vmin=0.0, vmax=1.0, aspect=1)
    imr = axr.imshow(dfr.values, cmap=cmap, vmin=0.0, vmax=1.0, aspect=1)

    axl.set_ylabel('Meiotic cell')
    axl.set_yticks([0, 1 ,2])
    axl.set_yticklabels(['Human', 'Mouse', 'Insect'])
    axl.xaxis.tick_top()
    axl.set_xticks([0, 1, 2])
    axl.set_xticklabels(['x endo.', 'x ecto.', 'x meso.'], rotation=30, ha='left')

    axr.set_yticks([0])
    axr.set_yticklabels([])
    axr.xaxis.tick_top()
    axr.set_xticks([0, 1])
    axr.set_xticklabels(['x MM meiotic', 'x DM meiotic'], rotation=30, ha='left')

    # Create colorbar
    cbar = plt.colorbar(iml, cax=cax, orientation='vertical', ticks=[0, 0.5, 1])
    cbar.ax.tick_params(labelsize='small')
    cbar.ax.set_ylabel('Similarity range', rotation=90, fontsize='small')
    cbar.ax.yaxis.set_label_position("left")

    plt.subplots_adjust(left=0.22, right=0.80, bottom=0.00, top=0.66, wspace=0.4, hspace=0.0)
    wIMGfile = 'images/img-jaccard-genes.pdf'
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()