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
mpl.rcParams['figure.titlesize'] = 'medium'
mpl.rcParams['axes.titlesize'] = 'small'
mpl.rcParams['axes.labelsize'] = 'small'
mpl.rcParams['xtick.labelsize'] = 'x-small'
mpl.rcParams['ytick.labelsize'] = 'x-small'
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['hatch.linewidth'] = 0.5
mpl.rcParams['hatch.color'] = '#969696'
import matplotlib.pyplot as plt
from data import *

#
# Label Scatter Plot
#
def labelScatterPlot(labels, ax, xs, ys, stdn, printLabel):
    df = pd.DataFrame({'xs': xs, 'ys': ys}, index=labels)
    stdx, stdy = df['xs'].std(), df['ys'].std()
    df = df.loc[ (df['xs'] >= stdx * stdn) | (df['xs'] <= -(stdx * stdn)) | (df['ys'] >= stdy * stdn) | (df['ys'] <= -(stdy * stdn)) ]

    df = df.sort_values('xs', ascending=False)
    textobjs = list()
    for label, data in df.iterrows():
        textobj = ax.text(data['xs'], data['ys'], label, fontsize=8, alpha=1, zorder=10, ha='center', va='center')
        textobjs.append(textobj)
    return textobjs
#
# Annotate Scatter Plot
#
def annotatePoint(label, ax, x, y, xpad, ypad):
    obj = ax.annotate(
        label,
        xy=(x, y), xytext=(x, y),
        xycoords='data', textcoords='data', ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
        #arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0.5', alpha = 0.2),
        fontsize=9, alpha=0.65
    )
    return obj

#
# Plot SVD
#
def plot_svd(celltype='spermatocyte', network='thr', threhsold=0.5, layer='DM'):
    threshold_str = str(threshold).replace('.', 'p')
    #
    data_cell = data_cells[celltype]
    modules = data_cell['modules-svd']['modules'][layer]

    #
    print('Plotting SVD for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
    rPCAFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rSFile = 'results/svd/{celltype:s}/svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-s.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    wIMGFile = 'images/svd/{celltype:s}/img-svd-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

    dfPCA = pd.read_csv(rPCAFile, index_col=0, encoding='utf-8')
    s = pd.read_csv(rSFile, squeeze=True, index_col=0, header=0, encoding='utf-8')
    # print(dfPCA)
    # print(s)
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(figsize=(8.5, 9), nrows=3, ncols=3)

    title = 'SVD on {celltype:s}-{network:s}-{threshold:s}-{layer:s} network'.format(celltype=celltype, network=network, threshold=str(threshold), layer=layer)
    fig.suptitle(title)

    # Plot EigenVals
    s_cumsum = s.cumsum()
    n_eigen_95 = s_cumsum[(s_cumsum < 0.95)].shape[0]

    n = 9
    ind = np.arange(n)
    height = s.iloc[:n].values
    width = 0.60
    xticklabels = (ind + 1)

    cmap = mpl.cm.get_cmap('hsv_r')
    norm = mpl.colors.Normalize(vmin=0, vmax=n)

    # S colors
    # tab20 = cm.get_cmap('tab20').colors
    # s_colors = tab20[0::2]
    # s_edgecolors = tab20[1::2]

    ax1.bar(ind, height, width, color='#1f77b4', edgecolor='#aec7e8', zorder=9, lw=1)
    ax1.set_xticks(ind)
    ax1.set_xticklabels(xticklabels)

    ax1.set_title('Explained variance ratio')
    ax1.annotate('95% with {:,d}\nsingular vectors'.format(n_eigen_95), xy=(0.97, 0.97), xycoords="axes fraction", ha='right', va='top', fontsize='x-small')
    ax1.set_xlabel('Components')
    ax1.set_ylabel('Variance')
    ax1.grid(axis='y')

    # Plot Projections
    for dim, ax in zip(range(1, 11), [ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]):
        print('- Dim: {:d}'.format(dim))
        #col = str(dim) + 'c'
        x = str(dim) + 'c'
        y = str(dim + 1) + 'c'
        xs = dfPCA[x].tolist()
        ys = dfPCA[y].tolist()
        #
        """
        if layer == 'DM':
            pca_colors = dfPCA['core'].map(lambda x: '#ff9896' if x == True else '#7f7f7f') # dfPCA['color'].tolist()
            pca_edgecolors = dfPCA['core'].map(lambda x: '#d62728' if x == True else '#c7c7c7') # dfPCA['color'].tolist()
        else:
            pca_colors = dfPCA['mammals'].map(lambda x: '#ff9896' if x == True else '#7f7f7f') # dfPCA['color'].tolist()
            pca_edgecolors = dfPCA['mammals'].map(lambda x: '#d62728' if x == True else '#c7c7c7') # dfPCA['color'].tolist()
        """
        pca_colors = '#7f7f7f'
        pca_edgecolors = '#c7c7c7'

        #if pca_colors == None:
        #   pca_colors = ys
        #   norm = mpl.colors.Normalize(vmin=min(ys), vmax=max(ys) )
        #   ax.scatter(xs,ys, c=pca_colors, cmap=cmapR2G2G, norm=norm, marker='o', edgecolor='black', lw=0.5, s=25)
        # else:
        ax.scatter(xs, ys, c=pca_colors, marker='o', edgecolor=pca_edgecolors, lw=0.5, s=20, zorder=5, rasterized=True)

        # Draw a X at the center
        ax.plot(0, 0, color='#2ca02c', marker='x', ms=16)

        # Draw lines at the center
        ax.axhline(y=0, c='black', lw=0.75, ls='-.', zorder=2)
        ax.axvline(x=0, c='black', lw=0.75, ls='-.', zorder=2)

        ax.set_title('Components {dim1} and {dim2}'.format(dim1=dim, dim2=(dim + 1)))
        ax.set_xlabel('Component {dim1:d}'.format(dim1=dim))
        ax.set_ylabel('Component {dim2:d}'.format(dim2=dim + 1))

        ax.grid()
        #ax.axis('equal')
        #ax.locator_params(axis='both', tight=True, nbins=6)

        #ax.set_ylim(-20,20)
        xlimmin, xlimmax = ax.get_xlim()
        ylimmin, ylimmax = ax.get_ylim()
        ylimdiff = abs(ylimmax) + abs(ylimmin)
        yperc = 0.035 * ylimdiff
        yspac = 0.8

        # Draw Component
        for module in modules:
            #
            xcomp = module['xy-coords']['x-comp']
            ycomp = module['xy-coords']['y-comp']
            if (dim == xcomp) and ((dim + 1) == ycomp):
                #
                name = "M{id:s}-{name:s}".format(id=str(module['id']), name=module['name'])
                facecolor = module.get('facecolor', 'black')
                edgecolor = module.get('edgecolor', 'none')
                hatch = module.get('hatch', None)
                x0, x1 = module['xy-coords']['x-values']
                y0, y1 = module['xy-coords']['y-values']
                # name loc
                name_loc  = module.get('name-loc', 'upper left')
                name_loc_upper_lower, name_loc_left_right = name_loc.split(' ')
                if name_loc_upper_lower == 'upper':
                    ytext = y1 + yperc
                    # add some space to y-lim-top
                    if abs(ylimmax - ytext) < 0.5:
                        ax.set_ylim((ylimmin - yspac, ylimmax + yspac))
                elif name_loc_upper_lower == 'lower':
                    ytext = y0 - yperc
                    # add some space to y-lim-bottom
                    if abs(ylimmin - ytext) < 0.5:
                        ax.set_ylim((ylimmin - yspac, ylimmax + yspac))
                if name_loc_left_right == 'left':
                    xtext = x0
                    ha = 'left'
                elif name_loc_left_right == 'right':
                    xtext = x1
                    ha = 'right'

                ax.fill([x0, x0, x1, x1], [y0, y1, y1, y0], facecolor=facecolor, edgecolor=edgecolor, lw=1, zorder=2, alpha=1, hatch=hatch)
                ann = ax.annotate(s=name, xy=(xtext, ytext), fontsize='xx-small', zorder=8, fontweight='bold', ha=ha, va='center')
                #           
            else:
                continue

    plt.subplots_adjust(left=0.07, right=0.98, bottom=0.06, top=0.93, wspace=0.30, hspace=0.35)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()


if __name__ == '__main__':

    celltype = 'enterocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    layer = 'HS'

    plot_svd(celltype, network, threshold, layer)
    #asd
    for celltype in ['spermatocyte', 'enterocyte']:
        for layer in ['HS', 'MM', 'DM']:
            plot_svd(celltype, network, threshold, layer)
