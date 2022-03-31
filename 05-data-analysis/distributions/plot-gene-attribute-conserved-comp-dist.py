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
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
#
mpl.rc('font', size=10)  # controls default text sizes
mpl.rc('axes', titlesize=12)  # fontsize of the axes title
mpl.rc('axes', labelsize=12)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=10)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=10)  # fontsize of the tick labels
mpl.rc('legend', fontsize=10)  # legend fontsize
mpl.rc('figure', titlesize=12)  # fontsize of the figure title
#
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
    parser.add_argument("--species", default="HS", type=str, help="Species. Defaults to 'HH'")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    biotype = args.biotype
    attribute = args.attribute
    log = args.log
    minTPM = args.minTPM
    species = layer = args.species

    print('Plot {celltype:s}-{biotype:s}-conserved comparison-{attribute:s} - {species:s}'.format(celltype=celltype, biotype=biotype, attribute=attribute, species=species))

    #
    # Conserved information
    #
    path_fpkm = '../../02-core_genes/results/FPKM/'
    df_HS = pd.read_csv(path_fpkm + 'HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
    df_MM = pd.read_csv(path_fpkm + 'MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
    df_DM = pd.read_csv(path_fpkm + 'DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')

    dict_string_gene_HS = df_HS['id_gene'].to_dict()
    dict_string_gene_MM = df_MM['id_gene'].to_dict()
    dict_string_gene_DM = df_DM['id_gene'].to_dict()

    print('Loading {celltype:s} meta genes'.format(celltype=celltype))
    path = '../../02-core_genes/results/'
    dfM = pd.read_csv(path + 'meta-genes/meta-{celltype:s}-genes.csv.gz'.format(celltype=celltype), index_col='id_eggnog', usecols=['id_eggnog', 'id_string_HS', 'id_string_MM', 'id_string_DM'])

    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    dfM['id_gene_HS'] = dfM['id_string_HS'].apply(lambda x: [dict_string_gene_HS[i] for i in x])
    dfM['id_gene_MM'] = dfM['id_string_MM'].apply(lambda x: [dict_string_gene_MM[i] for i in x])
    dfM['id_gene_DM'] = dfM['id_string_DM'].apply(lambda x: [dict_string_gene_DM[i] for i in x])

    dfM = dfM[['id_gene_HS', 'id_gene_MM', 'id_gene_DM']]
    # Only keep meta genes with homologs in all three species
    dfM = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 3]

    dict_conserved = {gene: True for gene in dfM['id_gene_' + species].explode().tolist()}
    #
    # Expression information
    #
    print('Loading {celltype:s} Files'.format(celltype=celltype))
    path = '../../02-core_genes/results/'
    dfS = pd.read_csv(path + 'FPKM/{species:s}/{species:s}-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype, species=species), index_col='id_gene')

    # Remove Duplicates
    dfS = dfS.loc[~dfS.index.duplicated(keep='first'), :]

    # minLogTPM
    if minTPM:
        dfS = dfS.loc[(dfS['TPM'] >= minTPM), :]

    if log is not None:
        dfS[attribute] = np.log2(dfS[attribute])
        #
        attribute_label = r'log$_2$({attribute:s})'.format(attribute=attribute)
    else:
        attribute_label = r'{attribute:s}'.format(attribute=attribute)

    if biotype is not None:
        print("Select nodes where biotype='{biotype:s}'".format(biotype=biotype))
        dfS = dfS.loc[dfS['biotype'] == biotype]

    # Adds Conserved information
    dfS['conserved'] = dfS.index.map(dict_conserved)

    # Separate conserved and not-conserved
    dfC = dfS.loc[dfS['conserved'] == True, :].copy()
    dfN = dfS.loc[dfS['conserved'] != True, :].copy()

    dict_celltype_translate = {
        'spermatocyte': 'meiotic',
        'spermatogonia': 'spermatogonia',
        'spermatid': 'spermatid',
        'enterocyte': 'intestine',
        'neuron': 'neuron',
        'muscle': 'muscle'
    }
    dict_species_color = {
        'HS': '#2ca02c',
        'MM': '#7f7f7f',
        'DM': '#ff7f0e',
    }
    celltype_str = dict_celltype_translate[celltype]
    title = r'All expressed {celltype:s} genes (TPM$\geq$1)'.format(celltype=celltype_str)
    facecolor = dict_species_color[species]
    facecolor = '#1f77b4'

    dfC.sort_values(attribute, ascending=False, inplace=True)
    dfN.sort_values(attribute, ascending=False, inplace=True)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    # title
    ax.set_title(title)

    # Plots
    max_value = math.ceil(max(dfC[attribute].max(), dfN[attribute].max()))
    bins = np.linspace(0, 10, 23, endpoint=True)
    hist_values_C = np.ones_like(dfC[attribute]) / len(dfC)
    hist_values_N = np.ones_like(dfN[attribute]) / len(dfN)
    (ac1, ac2, dc) = ax.hist(dfC[attribute], bins=bins, weights=hist_values_C, edgecolor='#d62728', facecolor=(214/255, 39/255, 40/255, 0.3), lw=1.5, zorder=5)
    (ad1, ad2, dn) = ax.hist(dfN[attribute], bins=bins, weights=hist_values_N, edgecolor='#1f77b4', facecolor=(31/255, 119/255, 180/255, 0.3), lw=1.5, zorder=4)

    # Insert Distribution Plot
    #axin = inset_axes(ax, width='40%', height='40%', loc='upper right', bbox_to_anchor=(.09, .19, .9, .8), bbox_transform=ax.transAxes)
    #pc, = axin.plot(np.arange(1, len(dfC) + 1), dfC[attribute], lw=1, marker='o', ms=3, markerfacecolor='none', markeredgecolor='#d62728', rasterized=True, zorder=5)
    #pn, = axin.plot(np.arange(1, len(dfN) + 1), dfN[attribute], lw=1, marker='o', ms=3, markerfacecolor='none', markeredgecolor=facecolor, rasterized=True, zorder=4)

    #axin.set_ylabel('Expression level'.format(label=attribute_label), fontsize=10)
    #axin.set_xlabel('Gene rank', fontsize=10)
    #ax.set_xscale('log')
    ax.set_ylabel('Probability')
    ax.set_xlabel(attribute_label)
    ax.set_xlim(-0.1, 10.1)

    # Legend
    ax.legend(
        handles=(dc, dn),
        labels=('Conserved', 'Non-conserved'),
        loc='upper right'
    )

    # Kolmogorov-Sminov test
    text = 'KS test\n'
    specie_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Fruit fly'}
    stat, pvalue = ks_2samp(dfC[attribute], dfN[attribute], alternative='two-sided', mode='auto')
    text += 'Conserved vs Not conserved: {stat:.2f} ({pvalue:.2e})'.format(stat=stat, pvalue=pvalue) + '\n'

    wMDfile = 'images/conserved-{attribute:s}/{celltype:s}/img-{celltype:s}-conserved-comp-{attribute:s}-{species:s}-dist.pdf.md'.format(celltype=celltype, attribute=attribute, species=species)
    ensurePathExists(wMDfile)
    with open(wMDfile, 'w') as file:
        file.write(text)

    # Grid
    #ax.grid(zorder=1)
    #axin.grid(zorder=1)

    plt.subplots_adjust(left=0.13, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
    wIMGfile = 'images/conserved-{attribute:s}/{celltype:s}/img-{celltype:s}-conserved-comp-{attribute:s}-{species:s}-dist.pdf'.format(celltype=celltype, attribute=attribute, species=species)
    ensurePathExists(wIMGfile)
    fig.savefig(wIMGfile)
    plt.close()
