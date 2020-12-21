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
#mpl.rcParams['axes.labelsize'] = 'small'
#mpl.rcParams['xtick.labelsize'] = 'x-small'
#mpl.rcParams['ytick.labelsize'] = 'x-small'
#mpl.rcParams['legend.fontsize'] = 'small'
#mpl.rcParams['hatch.linewidth'] = 0.5
#mpl.rcParams['hatch.color'] = '#969696'
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import shapely.geometry as sg
import descartes
from cycler import cycler
#
from data_spermatocyte_pca_modules_dm import spermatocyte_pca_modules_dm
from data_spermatocyte_pca_modules_mm import spermatocyte_pca_modules_mm
from data_spermatocyte_pca_modules_hs import spermatocyte_pca_modules_hs
#
from data_enterocyte_pca_modules_dm import enterocyte_pca_modules_dm
from data_enterocyte_pca_modules_mm import enterocyte_pca_modules_mm
from data_enterocyte_pca_modules_hs import enterocyte_pca_modules_hs


def plot_pca(celltype='spermatocyte', network='thr', threshold=0.5, layer='DM', modules=[]):
    """ Plot PCA """
    threshold_str = str(threshold).replace('.', 'p')
    #
    print('Plotting PCA for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
    rPCAFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rDiAnFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dian.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rEntFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-entropy.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rSFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-s.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

    df_pca = pd.read_csv(rPCAFile, index_col=0)
    df_dian = pd.read_csv(rDiAnFile, index_col=0)
    df_ent = pd.read_csv(rEntFile, index_col=0)
    s = pd.read_csv(rSFile, squeeze=True, index_col=0, header=0, encoding='utf-8')

    df_cp = df_ent.loc[df_ent['cut-rank'].notnull(), :].sort_values(['dim', 'cut-rank'])
    #
    cyc = (cycler(edgecolor=['#1f77b4', '#ff7f0e', '#2ca02c']) + cycler(linestyle=['solid', 'dashed', 'dotted']))()

    # Plot Variance
    s_cumsum = s.cumsum()
    n_eigen_95 = s_cumsum[(s_cumsum < 0.95)].shape[0]

    n = 9
    ind = np.arange(n)
    height = s.iloc[:n].values
    width = 0.60
    xticklabels = (ind + 1)

    fig, ax = plt.subplots(figsize=(3, 3))
    ax.bar(ind, height, width, color='#636363', edgecolor='#969696', zorder=9, lw=1)
    ax.set_xticks(ind)
    ax.set_xticklabels(xticklabels)

    species_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    title = '{species:s}'.format(species=species_name[layer])
    ax.set_title(title)

    ax.annotate('95% with {:,d}\nsingular vectors'.format(n_eigen_95), xy=(0.97, 0.97), xycoords="axes fraction", ha='right', va='top', fontsize='small')
    ax.set_xlabel('Components')
    ax.set_ylabel('Variance')
    ax.grid(axis='y')
    plt.subplots_adjust(left=0.21, right=0.97, bottom=0.17, top=0.89)
    #plt.tight_layout()
    wIMGFile = 'images/pca-entropy/{celltype:s}/{layer:s}/img-pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-variance.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()

    # Plot Projections
    for dim in range(1, 9):
        print('- Dim: {dim:d}'.format(dim=dim))
        # col = str(dim) + 'c'
        x = str(dim) + 'c'
        y = str(dim + 1) + 'c'
        xs = df_pca[x].tolist()
        ys = df_pca[y].tolist()
        #
        facecolors = '#c7c7c7'
        edgecolors = 'black'  # '#c7c7c7'

        fig, ax = plt.subplots(figsize=(3, 3))
        ax.scatter(xs, ys, c=facecolors, marker='o', edgecolor=edgecolors, lw=0.2, s=10, zorder=5, rasterized=True)

        # Draw a X at the center
        ax.plot(0, 0, color='#2ca02c', marker='x', ms=16)

        # Draw lines at the center
        ax.axhline(y=0, c='black', lw=0.75, ls='-.', zorder=2)
        ax.axvline(x=0, c='black', lw=0.75, ls='-.', zorder=2)

        ax.set_title('Components {dim1} and {dim2}'.format(dim1=dim, dim2=(dim + 1)))
        ax.set_xlabel('Component {dim1:d}'.format(dim1=dim))
        ax.set_ylabel('Component {dim2:d}'.format(dim2=dim + 1))

        ax.grid()

        xlimmin, xlimmax = ax.get_xlim()
        ylimmin, ylimmax = ax.get_ylim()
        ylimdiff = abs(ylimmax) + abs(ylimmin)
        yperc = 0.035 * ylimdiff
        yspac = 0.8

        # Radius Circles
        df_cp_tmp = df_cp.loc[(df_cp['dim'] == dim), :]
        sg_circles = {}
        for idx, cut_rank, radius in df_cp_tmp[['cut-rank', 'radius-start']].to_records():
            # Shapely Circle
            sg_circle = sg.Point(0, 0).buffer(radius)
            sg_circles[cut_rank] = sg_circle
            # Mpl Circle
            props = next(cyc)
            mpl_circle = mp.Circle((0, 0), radius=radius, facecolor='none', zorder=6, **props)
            ax.add_patch(mpl_circle)

        # Draw Component
        for module in modules:
            #
            xc = module['dim-coords']['xdim']
            yc = module['dim-coords']['ydim']
            if (dim == xc) and ((dim + 1) == yc):
                #
                mid = module['id']
                mname = module['name']

                # only rename a DM-M12
                if mname in dict_replace.keys():
                    mname = dict_replace.get(mname)
                #
                cx = '{xc:d}c'.format(xc=xc)  # label 1 component
                cy = '{yc:d}c'.format(yc=yc)  # label 2 component
                cxy = '{xc:d}c-{yc:d}c-dist'.format(xc=xc, yc=yc)  # label-1c-2c-dist
                #
                facecolor = module.get('facecolor', 'black')
                edgecolor = module.get('edgecolor', 'none')
                hatch = module.get('hatch', None)
                cxl, cxh = module['dim-coords']['xvals']
                cyl, cyh = module['dim-coords']['yvals']
                cut_rank = module['dim-coords']['radius-rank']
                sg_circle = sg_circles[cut_rank]
                # Radius of the circle
                cut_radius = df_ent.loc[((df_ent['dim'] == xc) & (df_ent['cut-rank'] == cut_rank)), 'radius-start'].squeeze()

                # Select points in module
                df_pca_tmp = df_pca.loc[
                    (
                        (df_pca[cx] >= cxl) & (df_pca[cx] <= cxh) & (df_pca[cy] >= cyl) & (df_pca[cy] <= cyh) & (df_dian[cxy] >= cut_radius)
                    ), ['gene', cx, cy]].copy()
                n = df_pca_tmp.shape[0]

                name = "M{mid:d}-{mname:s} (n={n:,d})".format(mid=mid, mname=mname, n=n)
                # name loc
                name_loc = module.get('name-loc', 'upper left')
                name_loc_upper_lower, name_loc_left_right = name_loc.split(' ')
                if name_loc_upper_lower == 'upper':
                    ytext = cyh + yperc
                    # add some space to y-lim-top
                    if abs(ylimmax - ytext) < 0.5:
                        ax.set_ylim((ylimmin - yspac, ylimmax + yspac))
                elif name_loc_upper_lower == 'lower':
                    ytext = cyl - yperc
                    # add some space to y-lim-bottom
                    if abs(ylimmin - ytext) < 0.5:
                        ax.set_ylim((ylimmin - yspac, ylimmax + yspac))
                if name_loc_left_right == 'left':
                    xtext = cxl
                    ha = 'left'
                elif name_loc_left_right == 'right':
                    xtext = cxh
                    ha = 'right'

                sg_box = sg.box(cxl, cyl, cxh, cyh)
                sg_diff = sg_box.difference(sg_circle)

                # ax.add_patch(descartes.PolygonPatch(sg_box, fc='b', ec='k', alpha=0.4))
                ax.add_patch(descartes.PolygonPatch(sg_diff, facecolor=facecolor, edgecolor=edgecolor, lw=1, zorder=2, alpha=1, hatch=hatch))

                # ax.fill([x0, x0, x1, x1], [y0, y1, y1, y0], facecolor=facecolor, edgecolor=edgecolor, lw=1, zorder=2, alpha=1, hatch=hatch)
                ax.annotate(text=name, xy=(xtext, ytext), fontsize='x-small', zorder=8, fontweight='bold', ha=ha, va='center')
            else:
                continue

        plt.subplots_adjust(left=0.21, right=0.97, bottom=0.17, top=0.89)
        #plt.tight_layout()
        wIMGFile = 'images/pca-entropy/{celltype:s}/{layer:s}/img-pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-C{dimx:d}x{dimy:d}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer, dimx=dim, dimy=(dim + 1))
        ensurePathExists(wIMGFile)
        plt.savefig(wIMGFile, dpi=300, bbox_inches=None, pad_inches=0.0)
        plt.close()


if __name__ == '__main__':

    # celltype = 'spermatocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    # layer = 'DM'

    data = {
        'spermatocyte': {
            'HS': spermatocyte_pca_modules_hs,
            'MM': spermatocyte_pca_modules_mm,
            'DM': spermatocyte_pca_modules_dm,
        },
        'enterocyte': {
            'HS': enterocyte_pca_modules_hs,
            'MM': enterocyte_pca_modules_mm,
            'DM': enterocyte_pca_modules_dm,
        }
    }

    dict_replace = {
        'Peptidyl-histidine dephosphorylation': 'Peptidyl-histidine dephospho.'
    }

    # modules = data[celltype][layer]
    # plot_pca(celltype, network, threshold, layer, modules)
    # asd
    for celltype in ['spermatocyte', 'enterocyte']:
        for layer in ['HS', 'MM', 'DM']:
            modules = data[celltype][layer]
            plot_pca(celltype, network, threshold, layer, modules)
