# coding=utf-8
# Author: Rion B Correia
# Date: Aug 06, 2019
#
# Description: Plots results of screened DM genes
#
# Instructions:
#
import numpy as np
import pandas as pd
import networkx as nx
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
from matplotlib import colors
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
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import product, chain


def value_to_color(x, cmap, norm):
    rgb = cmap(norm(x))[:3]
    return colors.rgb2hex(rgb)


def calc_control_mean_std_fert_rate(x):
    fertrate = x['hatched'] / x['eggs']
    return pd.Series({'mean fert-rate': fertrate.mean(), 'std fert-rate': fertrate.std()})


def get_network_layer(G, layer=''):
    return G.subgraph([n for n, d in G.nodes(data=True) if (d.get('layer') == layer)]).copy()


if __name__ == '__main__':

    # Load genes
    # df = pd.read_csv('../02-core_genes/results/pipeline-core/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])

    # Load Screened data
    dfs = pd.read_csv('data/core_DM_screened_2022-02-21.csv', index_col=0)

    # Load Control data
    dfc = pd.read_csv('data/screened_DM_controls.csv', index_col=0)
    dfc = dfc.groupby(dfc.index).apply(calc_control_mean_std_fert_rate)

    # Load FPKM data
    # dfFPKM = pd.read_csv('../02-core_genes/results/FPKM/DM/DM-FPKM-spermatocyte.csv.gz', index_col=0, usecols=['id_gene', 'FPKM'])

    dfs_only = dfs.loc[~dfs['FT1 eggs'].isnull(), :].copy()

    #status_cats = ['Screened', 'To be crossed', 'Pending', 'Reorder']
    #dfs['Status'] = pd.Categorical(dfs['Status'], categories=status_cats, ordered=True)
    #df['Status'] = dfs['Status']

    cols = ['gene', 'FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched']
    df = dfs_only[cols].copy()

    # Only plot screened genes
    #df = df.loc[df['Status'] == 'Screened', :]

    # Calculations
    df['total-eggs'] = 0
    df['total-hatched'] = 0
    for ft in range(1, 5):

        col_eggs = 'FT{:d} eggs'.format(ft)
        col_hatched = 'FT{:d} hatched'.format(ft)
        col_fertate = 'FT{:d} fert-rate'.format(ft)
        df[col_fertate] = df[col_hatched] / df[col_eggs]
        df['total-eggs'] += df[col_eggs]
        df['total-hatched'] += df[col_hatched]

    # Mean/SD
    df['mean fert-rate'] = df[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].mean(axis=1)
    df['std fert-rate'] = df[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].std(axis=1)

    #print(dfs.head())
    df['RNAi'] = dfs['Previous ref to RNAi working?']
    df['our-DM-code'] = dfs['Our DM pheno code']
    df['ext-DM-code'] = dfs['Others DM pheno code']
    df['ext-MM-code'] = dfs['MM pheno code']
    df['ext-HS-code'] = dfs['HS pheno code']

    # FPKM
    # df['FPKM'] = dfFPKM['FPKM']
    # df['logFPKM'] = df['FPKM'].apply(lambda x: np.log2(x + 1))

    # maxfpkm, minfpkm = df['logFPKM'].max(), df['logFPKM'].min()

    codemap = {

    }
    code_label = {
        'A': 'Meiotic',
        'B': 'Post-meiotic',
        'C': 'Gametes',
        'D': 'Pre-meiotic',
        #'E': 'General impairment of spermatogenesis',
        'F': 'Undetectable',
        #'G': 'Unspecified ',
        #'H': 'Non-germ cell autonomous'
    }
    code_facecolor = {
        'A': '#d62728',
        'B': '#e377c2',
        'C': '#9467bd',
        'D': '#2ca02c',
        #'E': '#9edae5',
        'F': '#f5d3a6',
        #'G': '#dadaeb',
        #'H': '#bdbdbd'
    }
    code_edgecolor = {
        'A': '#ff9896',
        'B': '#f7b6d2',
        'C': '#c5b0d5',
        'D': '#98df8a',
        #'E': '#9edae5',
        'F': '#fae7cf',
        #'G': '#dadaeb',
        #'H': '#bdbdbd'
    }

    # print(df.head())
    # print(df.tail())
    # print("logFPKM: {:.2f}/{:.2f}".format(minfpkm, maxfpkm))

    #
    # Plot all data
    #
    df = df.sort_values(['mean fert-rate', 'our-DM-code'], ascending=[True, True]).reset_index()

    print("Plotting")    
    fig = plt.figure(figsize=(7.5, 2.5))
    # fig.suptitle('Core metazoan meiotic genes'.format(page, number_of_pages))

    gs = gridspec.GridSpec(nrows=12, ncols=1)
    ax_fert = plt.subplot(gs[1:11, :1])

    adjustable = 'datalim'
    aspect = 'auto'
    ax_fert.set(adjustable=adjustable, aspect=aspect, anchor='NE')

    rotation = 50
    s = 5
    marker = '|'
    n = len(df)
    xticks = list(np.arange(0, n, 100))
    xticklabels = xticks #[::-1]

    #
    # DM Expression Values
    #
    for idx, row in df.iterrows():
        if row['mean fert-rate'] == 0:
            markerfacecolor = code_facecolor[row['our-DM-code']]
            markeredgecolor = code_edgecolor[row['our-DM-code']]
        elif row['mean fert-rate'] < 0.75:
            markerfacecolor = code_edgecolor[row['our-DM-code']]
            markeredgecolor = code_facecolor[row['our-DM-code']]
        elif row['mean fert-rate'] == 1.0:
            markerfacecolor = '#1f77b4'
            markeredgecolor = '#aec7e8'
        else:
            markerfacecolor = '#aec7e8'
            markeredgecolor = '#1f77b4'
        eb = ax_fert.errorbar(idx, row['mean fert-rate'], yerr=row['std fert-rate'], lw=0,
                              ecolor=markeredgecolor, elinewidth=0.6, capsize=0.6,
                              marker='.', markersize=1.5,
                              markeredgecolor=markeredgecolor, markeredgewidth=0.0,
                              markerfacecolor=markerfacecolor, markerfacecoloralt=None, zorder=5)

    """
    eb = ax_fert.errorbar(range(0, len(df)), df['mean fert-rate'], yerr=df['std fert-rate'], lw=0,
                          ecolor='#3182bd', elinewidth=0.6, capsize=0.6,
                          marker='.', markersize=1.5,
                          markeredgecolor='#3182bd', markeredgewidth=0.0,
                          markerfacecolor='#6baed6', markerfacecoloralt=None, zorder=5)
    """
    ax_fert.axhline(0.75, color='#d62728', lw=1, zorder=6)
    ax_fert.set_ylabel('Male fertility rate')
    ax_fert.set_xlabel('Candidate genes rank')
    ax_fert.set_yticks(np.linspace(0, 1, 5))
    # ax_fert.set_xticklabels([], fontsize='small', rotation=0)
    ax_fert.set_xticks(xticks)
    ax_fert.set_xticklabels(xticklabels)
    ax_fert.set_ylim(-0.04, 1.04)
    ax_fert.set_xlim(-1, len(df))
    ax_fert.grid(which='major', axis='y', linewidth=0.5)
    ax_fert.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax_fert.tick_params(axis='both')
    #ax_fert.invert_yaxis()

    # Red background
    ax_fert.fill_betweenx([-0.1, 1.1], [0, 0], [250, 250], color='gray', alpha=0.1, zorder=1)

    plt.subplots_adjust(left=0.10, right=0.98, bottom=0.13, top=0.999, wspace=0.0, hspace=0.0)
    fig.savefig('images/img-core_DM_screened-horizontal-meanfert.pdf')
    plt.close()
