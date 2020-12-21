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
import swifter


#  Separating by At Least One Match
def select_by_at_least_one_match(ilist, keeplist):
    # Only keep genes that are found in any of our gene list
    return [i for i in ilist if i in keeplist]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
    parser.add_argument("--biotype", default='protein_coding', type=str, choices=['protein_coding'], help="Filter nodes by biotype (e.g., protein-coding)")
    parser.add_argument("--attribute", default='TPM', type=str, help="Which attribute to plot. Defaults to 'TPM'.")
    parser.add_argument("--log", default=True, type=bool, help="Transforms attribute into log2(attribute).")
    parser.add_argument("--minTPM", default=1, type=int, help="minLogTPM = math.log2(x). Defaults to 1.")
    parser.add_argument("--title", default="Conserved meiotic genes")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    biotype = args.biotype
    attribute = args.attribute
    log = args.log
    minTPM = args.minTPM
    title = args.title

    dict_celltype_translate = {
        'spermatocyte': 'meiotic',
        'spermatogonia': 'spermatogonia',
        'spermatid': 'spermatid',
        'enterocyte': 'intestine',
        'neuron': 'neuron',
        'muscle': 'muscle'
    }
    celltype_str = dict_celltype_translate[celltype]

    print('Plot {celltype:s}-{biotype:s}-{attribute:s}'.format(celltype=celltype, biotype=biotype, attribute=attribute))

    print('Loading {celltype:s} Files'.format(celltype=celltype))
    path = '../../02-core_genes/results/'
    df_HS = pd.read_csv(path + 'FPKM/HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')
    df_MM = pd.read_csv(path + 'FPKM/MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')
    df_DM = pd.read_csv(path + 'FPKM/DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')
    dfM = pd.read_csv(path + 'meta-genes/meta-{celltype:s}-genes.csv.gz'.format(celltype=celltype), index_col='id_eggnog')

    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    # Remove Duplicates
    df_HS = df_HS.loc[~df_HS.index.duplicated(keep='first'), :]
    df_MM = df_MM.loc[~df_MM.index.duplicated(keep='first'), :]
    df_DM = df_DM.loc[~df_DM.index.duplicated(keep='first'), :]

    # minLogTPM
    if minTPM:
        df_HS = df_HS.loc[(df_HS['TPM'] >= minTPM), :]
        df_MM = df_MM.loc[(df_MM['TPM'] >= minTPM), :]
        df_DM = df_DM.loc[(df_DM['TPM'] >= minTPM), :]

    specie_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    specie_color = {'HS': '#2ca02c', 'MM': '#7f7f7f', 'DM': '#ff7f0e'}

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
        #
        title += ' (pc)'

    # Sort
    df_HS.sort_values(attribute, ascending=False, inplace=True)
    df_MM.sort_values(attribute, ascending=False, inplace=True)
    df_DM.sort_values(attribute, ascending=False, inplace=True)

    # Remove meta-genes not in selected species genes (df_HS) ..
    dfM['id_string_HS'] = dfM['id_string_HS'].swifter.apply(select_by_at_least_one_match, args=(df_HS['id_string'].to_list(), ))
    dfM['id_string_MM'] = dfM['id_string_MM'].swifter.apply(select_by_at_least_one_match, args=(df_MM['id_string'].to_list(), ))
    dfM['id_string_DM'] = dfM['id_string_DM'].swifter.apply(select_by_at_least_one_match, args=(df_DM['id_string'].to_list(), ))

    # Pairwise loop
    for (specie_i, df_i), (specie_j, df_j) in combinations([('HS', df_HS), ('MM', df_MM), ('DM', df_DM)], 2):
        #
        pair = specie_i + 'x' + specie_j
        print('Plotting: {pair:s}'.format(pair=pair))
        #
        name_i = specie_name[specie_i]
        name_j = specie_name[specie_j]
        color_i = specie_color[specie_i]
        color_j = specie_color[specie_j]

        title = r'Conserved {celltype:s} genes {name_i:s} x {name_j:s}'.format(celltype=celltype_str, name_i=name_i, name_j=name_j)
        #
        axin_bbox_to_anchor = (.09, .19, .9, .8)
        axin_loc = 'upper right'
        legend_loc = 'lower left'

        # tmp dfM
        dfMt = dfM[['id_string_' + specie_i, 'id_string_' + specie_j]].copy()
        #
        # Index only rows with at least one id_string in each specie
        dfMt = dfMt.loc[dfMt.applymap(len).applymap(bool).sum(axis='columns') == 2]

        set_id_string_i = set(dfMt['id_string_' + specie_i].explode().to_list())
        set_id_string_j = set(dfMt['id_string_' + specie_j].explode().to_list())
        #
        df_i = df_i.loc[df_i['id_string'].isin(set_id_string_i), :]
        df_j = df_j.loc[df_j['id_string'].isin(set_id_string_j), :]

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

        # Title
        ax.set_title(title)

        # Plots
        pi, = ax.plot(np.arange(1, len(df_i) + 1), df_i[attribute], lw=0, marker='o', ms=4, color=color_i, rasterized=True, zorder=5)
        pj, = ax.plot(np.arange(1, len(df_j) + 1), df_j[attribute], lw=0, marker='o', ms=4, color=color_j, rasterized=True, zorder=4)

        max_value = math.ceil(max(df_i[attribute].max(), df_j[attribute].max()))
        bins = np.linspace(0, max_value, 23, endpoint=True)

        # Insert Distribution Plot
        hist_values_i = np.ones_like(df_i[attribute]) / len(df_i)
        hist_values_j = np.ones_like(df_j[attribute]) / len(df_j)
        axin = inset_axes(ax, width='40%', height='40%', loc=axin_loc, bbox_to_anchor=axin_bbox_to_anchor, bbox_transform=ax.transAxes)
        axin.hist(df_i[attribute], bins=bins, weights=hist_values_i, edgecolor=color_i, facecolor=(0, 0, 0, 0), lw=1, zorder=5)
        axin.hist(df_j[attribute], bins=bins, weights=hist_values_j, edgecolor=color_j, facecolor=(0, 0, 0, 0), lw=1, zorder=4)

        ax.set_ylabel('Expression level [{label:s}]'.format(label=attribute_label))
        ax.set_xlabel('Gene rank')
        ax.set_xscale('log')
        axin.set_ylabel('Probability', fontsize='small')
        axin.set_xlabel(attribute_label, fontsize='small')

        # Legend
        ax.legend(
            handles=(pi, pj),
            labels=(name_i, name_j),
            loc=legend_loc
        )

        # Kolmogorov-Sminov test
        text = 'KS test\n'
        stat, pvalue = ks_2samp(df_i[attribute], df_j[attribute], alternative='two-sided', mode='auto')
        text += '{i:s}-{j:}: {stat:.2f} ({pvalue:.2e})'.format(i=name_i, j=name_j, stat=stat, pvalue=pvalue) + '\n'
        #
        wMDfile = 'images/{attribute:s}/{celltype:s}/img-{celltype:s}-{attribute:s}-conserved-{pair:s}-dist.pdf.md'.format(celltype=celltype, attribute=attribute, pair=pair)
        ensurePathExists(wMDfile)
        with open(wMDfile, 'w') as file:
            file.write(text)
        #
        ax.annotate(s='KS={stat:.2f}'.format(stat=stat), xy=(0.9, 0.08), xytext=(0.9, 0.08), xycoords='axes fraction', ha='right', zorder=10)

        # Grid
        ax.grid(zorder=1)
        axin.grid(zorder=1)

        plt.subplots_adjust(left=0.13, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
        wIMGfile = 'images/{attribute:s}/{celltype:s}/img-{celltype:s}-{attribute:s}-conserved-{pair:s}-dist.pdf'.format(celltype=celltype, attribute=attribute, pair=pair)
        ensurePathExists(wIMGfile)
        fig.savefig(wIMGfile)
        plt.close()
