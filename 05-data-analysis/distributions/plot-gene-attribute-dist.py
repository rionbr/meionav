# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots full conserved network TPM (node attribute) distribution for all three species.
#
# Instructions:
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
import argparse
from utils import ensurePathExists
from scipy.stats import ks_2samp
from itertools import combinations


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
    parser.add_argument("--biotype", default='protein_coding', type=str, choices=['protein_coding'], help="Filter nodes by biotype (e.g., protein-coding)")
    parser.add_argument("--attribute", default='TPM', type=str, help="Which attribute to plot. Defaults to 'TPM'.")
    parser.add_argument("--log", default=True, type=bool, help="Transforms attribute into log2(attribute).")
    parser.add_argument("--minTPM", default=1, type=int, help="minLogTPM = math.log2(x). Defaults to 1.")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    biotype = args.biotype
    attribute = args.attribute
    log = args.log
    minTPM = args.minTPM

    print('Plot {celltype:s}-{biotype:s}-{attribute:s}'.format(celltype=celltype, biotype=biotype, attribute=attribute))

    print('Loading {celltype:s} Files'.format(celltype=celltype))
    path = '../../02-core_genes/results/'
    df_HS = pd.read_csv(path + 'FPKM/HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')
    df_MM = pd.read_csv(path + 'FPKM/MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')
    df_DM = pd.read_csv(path + 'FPKM/DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')

    # Remove Duplicates
    df_HS = df_HS.loc[~df_HS.index.duplicated(keep='first'), :]
    df_MM = df_MM.loc[~df_MM.index.duplicated(keep='first'), :]
    df_DM = df_DM.loc[~df_DM.index.duplicated(keep='first'), :]

    # minLogTPM
    if minTPM:
        df_HS = df_HS.loc[(df_HS['TPM'] >= minTPM), :]
        df_MM = df_MM.loc[(df_MM['TPM'] >= minTPM), :]
        df_DM = df_DM.loc[(df_DM['TPM'] >= minTPM), :]

    if log is not None:
        df_HS[attribute] = np.log2(df_HS[attribute])
        df_MM[attribute] = np.log2(df_MM[attribute])
        df_DM[attribute] = np.log2(df_DM[attribute])
        #
        attribute_label = r'log$_2$({attribute:s})'.format(attribute=attribute)
    else:
        attribute_label = r'{attribute:s}'.format(attribute=attribute)

    if biotype is not None:
        print("Select nodes where biotype='{biotype:s}'".format(biotype=biotype))
        df_HS = df_HS.loc[df_HS['biotype'] == biotype]
        df_MM = df_MM.loc[df_MM['biotype'] == biotype]
        df_DM = df_DM.loc[df_DM['biotype'] == biotype]

    dict_celltype_translate = {
        'spermatocyte': 'meiotic',
        'spermatogonia': 'spermatogonia',
        'spermatid': 'spermatid',
        'enterocyte': 'intestine',
        'neuron': 'neuron',
        'muscle': 'muscle'
    }
    celltype_str = dict_celltype_translate[celltype]
    title = r'All expressed {celltype:s} genes (TPM$\geq$1)'.format(celltype=celltype_str)

    df_HS.sort_values(attribute, ascending=False, inplace=True)
    df_MM.sort_values(attribute, ascending=False, inplace=True)
    df_DM.sort_values(attribute, ascending=False, inplace=True)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # Title
    ax.set_title(title)

    # Plots
    phs, = ax.plot(np.arange(1, len(df_HS) + 1), df_HS[attribute], lw=0, marker='o', ms=4, color='#2ca02c', rasterized=True, zorder=5)
    pmm, = ax.plot(np.arange(1, len(df_MM) + 1), df_MM[attribute], lw=0, marker='o', ms=4, color='#7f7f7f', rasterized=True, zorder=4)
    pdm, = ax.plot(np.arange(1, len(df_DM) + 1), df_DM[attribute], lw=0, marker='o', ms=4, color='#ff7f0e', rasterized=True, zorder=3)

    max_value = math.ceil(max(df_HS[attribute].max(), df_MM[attribute].max(), df_DM[attribute].max()))
    bins = np.linspace(0, max_value, 23, endpoint=True)

    # Insert Distribution Plot
    hist_values_HS = np.ones_like(df_HS[attribute]) / len(df_HS)
    hist_values_MM = np.ones_like(df_MM[attribute]) / len(df_MM)
    hist_values_DM = np.ones_like(df_DM[attribute]) / len(df_DM)
    axin = inset_axes(ax, width='40%', height='40%', loc='upper right', bbox_to_anchor=(.09, .19, .9, .8), bbox_transform=ax.transAxes)
    axin.hist(df_HS[attribute], bins=bins, weights=hist_values_HS, edgecolor='#2ca02c', facecolor=(0, 0, 0, 0), lw=1, zorder=5)
    axin.hist(df_MM[attribute], bins=bins, weights=hist_values_MM, edgecolor='#7f7f7f', facecolor=(0, 0, 0, 0), lw=1, zorder=4)
    axin.hist(df_DM[attribute], bins=bins, weights=hist_values_DM, edgecolor='#ff7f0e', facecolor=(0, 0, 0, 0), lw=1, zorder=3)

    ax.set_ylabel('Expression level [{label:s}]'.format(label=attribute_label))
    ax.set_xlabel('Gene rank')
    ax.set_xscale('log')
    axin.set_ylabel('Probability', fontsize='small')
    axin.set_xlabel(attribute_label, fontsize='small')

    # Legend
    ax.legend(
        handles=(phs, pmm, pdm),
        labels=('Human', 'Mouse', 'Insect'),
        loc='lower left'
    )

    # Kolmogorov-Sminov test
    text = 'KS test\n'
    specie_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    for (specie_i, vals_i), (specie_j, vals_j) in combinations([('HS', df_HS[attribute]), ('MM', df_MM[attribute]), ('DM', df_DM[attribute])], 2):
        stat, pvalue = ks_2samp(vals_i, vals_j, alternative='two-sided', mode='auto')
        text += '{i:s}-{j:}: {stat:.2f} ({pvalue:.2e})'.format(i=specie_name[specie_i], j=specie_name[specie_j], stat=stat, pvalue=pvalue) + '\n'

    wMDfile = 'images/{attribute:s}/{celltype:s}/img-{celltype:s}-{attribute:s}-dist.pdf.md'.format(celltype=celltype, attribute=attribute)
    ensurePathExists(wMDfile)
    with open(wMDfile, 'w') as file:
        file.write(text)

    # Grid
    ax.grid(zorder=1)
    axin.grid(zorder=1)

    plt.subplots_adjust(left=0.13, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
    wIMGfile = 'images/{attribute:s}/{celltype:s}/img-{celltype:s}-{attribute:s}-dist.pdf'.format(celltype=celltype, attribute=attribute)
    ensurePathExists(wIMGfile)
    fig.savefig(wIMGfile)
    plt.close()
