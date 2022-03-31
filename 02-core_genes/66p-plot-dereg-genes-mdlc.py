# coding=utf-8
# Author: Rion B Correia
# Date: Sept 15, 2021
#
# Description: Plot core deregulated genes in mdlc
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
#
mpl.rc('font', size=10)  # controls default text sizes
mpl.rc('axes', titlesize=12)  # fontsize of the axes title
mpl.rc('axes', labelsize=12)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=12)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=12)  # fontsize of the tick labels
mpl.rc('legend', fontsize=10)  # legend fontsize
mpl.rc('figure', titlesize=12)  # fontsize of the figure title
#
from matplotlib.colors import Normalize
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from utils import ensurePathExists


if __name__ == '__main__':

    dfs = pd.read_csv('../03-screened-data/data/core_DM_screened_2021-06-07.csv', index_col=0)
    # Only those that gave a phenotype
    dfs = dfs.loc[dfs['Recorded cellular phenotype'].notnull(), :]

    # mdlc vs controls
    df = pd.read_csv('../01-diff-gene-exp/results/mdlc/DM-DGE-mdlc_vs_control.csv', index_col=0)

    # Only up/down regulated
    df = df.loc[(df['FDR'] <= 0.05) & (df['logCPM'] > 1) & (df['logFC'].abs() > 1), :]

    # Only CORE
    #df = df.loc[df.index.isin(dfc.index), :]
    # Only PHENOTYPE
    df = df.loc[df.index.isin(dfs.index), :].copy()

    ccols = ['mC_1', 'mC_2', 'mC_3']
    mcols = ['mdlc_1', 'mdlc_2', 'mdlc_3']
    cols = ccols + mcols
    df[cols] = df[cols].apply(zscore, axis='columns', raw=True)

    #vmax, vmin = df[cols].max().max(), df[cols].min().min()

    # Sort
    df.sort_values('logFC', inplace=True)
    #
    df_u = df.loc[(df['logFC'] < 0), :]
    df_d = df.loc[(df['logFC'] > 0), :]
    #
    df_cu = df_u[ccols]
    df_mu = df_u[mcols]
    df_cd = df_d[ccols]
    df_md = df_d[mcols]

    #
    # Plot Entry
    #
    fig = plt.figure(figsize=(4.0, 9.4))
    #
    gs = GridSpec(nrows=46, ncols=2)
    ax_cu = plt.subplot(gs[0:30, 0])
    ax_cd = plt.subplot(gs[31:46, 0])
    ax_mu = plt.subplot(gs[0:30, 1])
    ax_md = plt.subplot(gs[31:46, 1])
    #
    dfl = [df_cu, df_mu, df_cd, df_md]
    axes = [ax_cu, ax_mu, ax_cd, ax_md]
    #
    # red = '#d62728'
    # blue = '#1f77b4'
    # green = '#2ca02c'
    cmap = 'coolwarm'
    vmin, vmax = -1.4, 1.4
    norm = Normalize(vmin=vmin, vmax=vmax)
    aspect_u = aspect_d = 'auto'

    # Control
    im = ax_cu.imshow(df_cu, cmap=cmap, norm=norm, aspect=aspect_u, interpolation='nearest')
    im = ax_cd.imshow(df_cd, cmap=cmap, norm=norm, aspect=aspect_d, interpolation='nearest')
    # mdlc
    im = ax_mu.imshow(df_mu, cmap=cmap, norm=norm, aspect=aspect_u, interpolation='nearest')
    im = ax_md.imshow(df_md, cmap=cmap, norm=norm, aspect=aspect_d, interpolation='nearest')

    ax_cu.set_title('Control (RNAi)',)
    ax_mu.set_title('       dRNF113 RNAi #1')

    ax_cu.set_yticks([])
    ax_cd.set_yticks([])
    ax_cu.set_xticks([])
    ax_mu.set_xticks([])

    ax_cd.set_xticks([0, 1, 2])
    ax_md.set_xticks([0, 1, 2])
    ax_cd.set_xticklabels(['R1', 'R2', 'R3'])
    ax_md.set_xticklabels(['R1', 'R2', 'R3'])

    ax_mu.set_yticks(range(len(df_mu)))
    ax_md.set_yticks(range(len(df_md)))
    ax_mu.set_yticklabels(df_u['gene'])
    ax_md.set_yticklabels(df_d['gene'])

    ax_mu.yaxis.tick_right()
    ax_md.yaxis.tick_right()
    # Colorbar
    cax = plt.axes([0.88, 0.03, 0.022, 0.12])
    cbd = plt.colorbar(im, orientation='vertical', cax=cax)
    cax.set_ylabel('z-score')
    cax.yaxis.set_label_position("left")

    for tick in ax_mu.get_yticklabels():
        if tick._text == 'Ubi-p63E':
            tick.set_fontweight('bold')
        tick.set_color('#1f77b4')

    for tick in ax_md.get_yticklabels():
        tick.set_color('#d62728')

    for ax, dft in zip(axes, dfl):
        n, m = dft.shape
        ax.set_xticks(np.arange(m + 1) - .5, minor=True)
        ax.set_yticks(np.arange(n + 1) - .5, minor=True)
        ax.tick_params(axis="both", which='minor', length=0)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=1)

    plt.subplots_adjust(left=0.04, right=0.60, bottom=0.03, top=0.95, wspace=0.20, hspace=0.10)

    wIMGfile = 'images/dereg-mdlc/img-dereg-mdlc.pdf'
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=300)
    plt.close()
