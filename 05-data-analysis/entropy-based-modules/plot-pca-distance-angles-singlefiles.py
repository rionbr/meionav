# coding=utf-8
# Author: Rion B Correia
# Date: June 17, 2020
#
# Description: Plots the entropy
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
import matplotlib.pyplot as plt
from cycler import cycler


def plot_distance_and_angles(celltype, network, threshold, layer, radius_window, radius_overlap, angle_window, angle_overlap):
    """ Plots Distance and Angles """
    threshold_str = str(threshold).replace('.', 'p')

    print('Plotting Distance & Angles for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
    rDiAnFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dian.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rEntrFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-entropy.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    #
    df_dian = pd.read_csv(rDiAnFile, index_col=0)
    df_ent = pd.read_csv(rEntrFile, index_col=0)
    df_cp = df_ent.loc[df_ent['cut-rank'].notnull(), :].sort_values(['dim', 'cut-rank'])

    #
    cyc = (cycler(color=['#1f77b4', '#ff7f0e', '#2ca02c']) + cycler(linestyle=['solid', 'dashed', 'dotted']))()
    #
    for dim in range(1, 10):
        print('- Dim: {:d}'.format(dim))
        df_ent_tmp = df_ent.loc[df_ent['dim'] == dim].copy()
        df_cp_tmp = df_cp.loc[(df_cp['dim'] == dim), :]
        #
        fig, ax = plt.subplots(figsize=(3.66, 3))
        axt = ax.twinx()
        #
        cx = str(dim) + 'c'
        cy = str(dim + 1) + 'c'
        dist_label = '{cx:s}-{cy:s}-dist'.format(cx=cx, cy=cy)
        angle_label = '{cx:s}-{cy:s}-angle'.format(cx=cx, cy=cy)
        facecolors = '#c7c7c7'
        edgecolors = 'black'
        df_dian = df_dian.sort_values([dist_label, angle_label], ascending=[True, True])
        xs = df_dian[dist_label]
        ys = df_dian[angle_label]
        #
        #
        ax.scatter(xs, ys, c=facecolors, marker='o', edgecolors=edgecolors, lw=0.2, s=10, zorder=4, rasterized=True)
        axt.plot(df_ent_tmp['radius-start'], df_ent_tmp['entropy-norm'], color='#d62728', zorder=6, marker='.', markersize=3, lw=0)
        axt.plot(df_ent_tmp['radius-start'], df_ent_tmp['entropy-smooth'], color='#ff9896', zorder=5)

        # Plot Cut Points
        for idx, cut_rank, radius in df_cp_tmp[['cut-rank', 'radius-start']].to_records():
            props = next(cyc)
            ax.axvline(x=radius, zorder=6, **props)
        #
        ax.set_title('Components {dim1} and {dim2}'.format(dim1=dim, dim2=(dim + 1)))
        ax.set_xlabel(r'radius ($\theta_w = {radius_window:.2f}, \theta_o = {radius_overlap:.2f}$)'.format(radius_window=radius_window, radius_overlap=radius_overlap))
        ax.set_ylabel(r'angle ($\varphi_w = {angle_window:d}, \varphi_o = {angle_overlap:d}$)'.format(angle_window=angle_window, angle_overlap=angle_overlap))
        yticks = [-180, -135, -90, -45, 0, 45, 90, 135, 180]
        #yticklabels = [r'$-\pi(180\degree)$', r'$-\frac{3\pi}{4}(135\degree)$', r'$-\frac{\pi}{2}(90\degree)$', r'$-\frac{\pi}{4}(45\degree)$', r'$0(0)$', r'$\frac{\pi}{4}(45\degree)$', r'$\frac{\pi}{2}(90\degree)$', r'$\frac{3\pi}{4}(135\degree)$', r'$\pi(180\degree)$']
        yticklabels = [r'$-\pi$', r'$-\frac{3\pi}{4}$', r'$-\frac{\pi}{2}$', r'$-\frac{\pi}{4}$', r'$0$', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{4}$', r'$\pi$']
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize='medium')
        #
        axt.set_ylabel('entropy (normed)')
        axt.set_ylim(0, 1)
        ax.grid()
        
        #
        plt.subplots_adjust(left=0.17, right=0.84, bottom=0.17, top=0.89)
        #plt.tight_layout()
        wIMGFile = 'images/pca-entropy/{celltype:s}/{layer:s}/img-entropy-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-C{dimx:d}x{dimy:d}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer, dimx=dim, dimy=(dim + 1))
        ensurePathExists(wIMGFile)
        plt.savefig(wIMGFile, dpi=300)
        plt.close()


if __name__ == '__main__':

    # celltype = 'spermatocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    # layer = 'HS'

    # from calc-pca-entropy.py
    radius_window = 1.0
    radius_overlap = 0.1
    angle_window = 30
    angle_overlap = 15

    #plot_distance_and_angles(celltype, network, threshold, layer, radius_window, radius_overlap, angle_window, angle_overlap)

    for celltype in ['spermatocyte', 'enterocyte']:
        for layer in ['HS', 'MM', 'DM']:
            plot_distance_and_angles(celltype, network, threshold, layer, radius_window, radius_overlap, angle_window, angle_overlap)
