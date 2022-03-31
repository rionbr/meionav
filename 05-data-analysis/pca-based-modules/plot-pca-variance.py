# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and computed its modules using Louvain & Infomap.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
#mpl.rcParams['figure.titlesize'] = 'medium'
#mpl.rcParams['axes.titlesize'] = 'small'
mpl.rcParams['axes.labelsize'] = 'small'
mpl.rcParams['xtick.labelsize'] = 'x-small'
mpl.rcParams['ytick.labelsize'] = 'x-small'
mpl.rcParams['legend.fontsize'] = 'small'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from data import *


def bar_labels(ax, rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        rect_x = rect.get_x()
        rect_y = rect.get_y()
        rect_height = rect.get_height()
        rect_width = rect.get_width()
        percent = rect_height
        pxs = [rect_x + rect_width + 0.075] * 2
        pys = [rect_y, rect_y + rect_height]
        ax.plot(pxs, pys, color='#d62728')
        ax.annotate('{percent:.1%}'.format(percent=percent),
                    xy=(rect_x + rect_width + 0.12, rect_y + rect_height / 2),
                    xytext=(5, 0),  # 3 points vertical offset
                    textcoords="offset points",
                    color='#d62728',
                    fontsize='small',
                    rotation=90,
                    ha='center', va='center')


if __name__ == '__main__':

    celltype = 'enterocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5

    threshold_str = str(threshold).replace('.', 'p')
    #
    wIMGFile = 'images/{celltype:s}/img-pca-variance-{celltype:s}-{network:s}-{threshold:s}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str)
    variances = []

    fig = plt.figure(figsize=(8.5, 2.4))
    gs = gridspec.GridSpec(ncols=18, nrows=1, figure=fig)
    ax_hs = fig.add_subplot(gs[0, 0:4])
    ax_mm = fig.add_subplot(gs[0, 5:9])
    ax_dm = fig.add_subplot(gs[0, 10:14])
    ax_bar = fig.add_subplot(gs[0, 15:18])

    species_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    var_ratios = {'HS': {'first9': 0, 'others': 0}, 'MM': {'first9': 0, 'others': 0}, 'DM': {'first9': 0, 'others': 0}}

    for layer, ax in zip(['HS', 'MM', 'DM'], [ax_hs, ax_mm, ax_dm]):
        print('Plotting PCA variance for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
        rSFile = '../../04-network/results/pca/{celltype:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-s.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

        s = pd.read_csv(rSFile, squeeze=True, index_col=0, header=0, encoding='utf-8')

        specie = species_name[layer]
        title = '{specie:s}'.format(specie=specie)
        ax.set_title(title)

        # Plot EigenVals
        s_cumsum = s.cumsum()
        n_eigen_95 = s_cumsum[(s_cumsum < 0.95)].shape[0]

        # Data for last axes
        var_ratios[layer]['first9'] = s[:9].sum()
        var_ratios[layer]['others'] = s[9:].sum()

        n = 9
        ind = np.arange(n)
        height = s.iloc[:n].values
        width = 0.60
        xticklabels = (ind + 1)

        ax.bar(ind, height, width, color='#1f77b4', edgecolor='#aec7e8', zorder=9, lw=1)
        ax.set_xticks(ind)
        ax.set_xticklabels(xticklabels)

        ax.annotate('95% with {:,d}\nsingular vectors'.format(n_eigen_95), xy=(0.97, 0.97), xycoords="axes fraction", ha='right', va='top', fontsize='x-small')
        ax.set_xlabel('Components')
        ax.set_ylabel('Variance')
        ax.grid(axis='y')

    ind = np.arange(3)
    width = 0.35
    first9 = np.array([var_ratios['HS']['first9'], var_ratios['MM']['first9'], var_ratios['DM']['first9']])
    others = np.array([var_ratios['HS']['others'], var_ratios['MM']['others'], var_ratios['DM']['others']])
    # first9 = first9 + others
    rects_first9 = ax_bar.bar(ind, first9, width, bottom=others, color='#1f77b4', edgecolor='#aec7e8', lw=1, zorder=9)
    rects_others = ax_bar.bar(ind, others, width, bottom=0, color='#bdbdbd', edgecolor='#d9d9d9', lw=1, zorder=9)

    bar_labels(ax=ax_bar, rects=rects_first9)

    ax_bar.set_xticks(ind)
    ax_bar.set_xticklabels(['Human', 'Mouse', 'Insect'], fontsize='small')
    ax_bar.set_ylabel('Variance')
    ax_bar.set_xlim(-0.4, 2.7)
    #ax_bar.grid(axis='y')

    plt.subplots_adjust(left=0.06, right=0.98, bottom=0.18, top=0.86, wspace=0.5, hspace=0.0)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()
