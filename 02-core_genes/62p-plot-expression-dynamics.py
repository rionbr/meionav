# coding=utf-8
# Author: Rion B Correia
# Date: Nov 05, 2020
#
# Description: Plots an overview of the FPKM across the three species
#
#
import numpy as np
import math
import pandas as pd
idx = pd.IndexSlice
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from utils import ensurePathExists
from sklearn import preprocessing


def load_species_fpkm_dfs(species='DM'):
    ldfs = []
    minTPM = 1
    #
    for celltype in ['spermatogonia', 'spermatocyte', 'spermatid']:
        print('Loading: {species:s} - {celltype:s}'.format(species=species, celltype=celltype))

        rCSVFile = 'results/FPKM/{species:s}/{species:s}-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype, species=species)
        df = pd.read_csv(rCSVFile, index_col=0)

        # Only TPM > 1
        df = df.loc[(df['TPM'] >= minTPM), :]
        # Drop Duplicated
        df = df.loc[~df.index.duplicated(), :]
        # Calc Log(FPKM)
        df['logFPKM'] = df['FPKM'].apply(lambda x: np.log2(x + 1))
        # Index to MultiIndex
        index = pd.MultiIndex.from_arrays([[celltype]*6, df.columns])
        df.columns = index

        ldfs.append(df)
    dfc = pd.concat(ldfs, axis='columns', verify_integrity=True).fillna(0)
    dfc.sort_values([('spermatogonia', 'logFPKM'), ('spermatocyte', 'logFPKM'), ('spermatid', 'logFPKM')], ascending=[False, False, False], inplace=True)

    return dfc


def plot_species_fpkm(df, species='DM'):
    fig, ax = plt.subplots(figsize=(2.5, 8))
    cax = plt.axes([0.12, 0.05, 0.04, 0.2])

    species_label = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    # red = '#d62728'
    # blue = '#1f77b4'
    # green = '#2ca02c'
    cmap = 'Reds'  # LinearSegmentedColormap.from_list(name='custom-cmap', colors=['white', '#d62728'])
    vmin = 0
    vmax = 10  # df.loc[:, idx[:, 'logFPKM']].max().max()
    norm = Normalize(vmin=vmin, vmax=vmax)

    dft = df.loc[:, idx[:, 'logFPKM']]
    im = ax.imshow(dft, cmap=cmap, norm=norm, aspect='auto', interpolation='nearest')

    #
    ax.yaxis.set_ticks_position('right')
    ax.set_title(species_label[species])
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['Gonia', 'Cyte', 'Tid'])

    plt.colorbar(im, cax=cax, extendfrac=0.1, extend='max')
    cax.yaxis.set_ticks_position('left')
    cax.set_ylabel('log(FPKM+1)')

    plt.subplots_adjust(left=0.29, right=0.77, bottom=0.05, top=0.96)
    # plt.tight_layout()
    wIMGfile = 'images/exp-dyn/img-exp-dyn-{species:s}.pdf'.format(species=species)
    plt.savefig(wIMGfile, dpi=300)


if __name__ == '__main__':

    df_dm = load_species_fpkm_dfs(species='DM')
    df_mm = load_species_fpkm_dfs(species='MM')
    df_hs = load_species_fpkm_dfs(species='HS')
    #
    plot_species_fpkm(df_dm, species='DM')
    plot_species_fpkm(df_mm, species='MM')
    plot_species_fpkm(df_hs, species='HS')
