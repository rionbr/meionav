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
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
from matplotlib import colors
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def value_to_color(x, cmap, norm):
    rgb = cmap(norm(x))[:3]
    return colors.rgb2hex(rgb)


def calc_control_mean_std_fert_rate(x):
    fertrate = x['hatched'] / x['eggs']
    return pd.Series({'mean fert-rate': fertrate.mean(), 'std fert-rate': fertrate.std()})


if __name__ == '__main__':

    rnai = 'Yes'  # None, 'Yes', 'No'
    # Load genes
    df = pd.read_csv('../02-core_genes/results/pipeline-core/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])

    # Load Screened data
    dfs = pd.read_csv('data/core_DM_screened_2020-10-21.csv', index_col=0)

    # Load Control data
    dfc = pd.read_csv('data/screened_DM_controls.csv', index_col=0)
    dfc = dfc.groupby(dfc.index).apply(calc_control_mean_std_fert_rate)

    dfs_only = dfs.loc[~dfs['FT1 eggs'].isnull(), :]

    status_cats = ['Screened', 'To be crossed', 'Pending', 'Reorder']
    dfs['Status'] = pd.Categorical(dfs['Status'], categories=status_cats, ordered=True)
    df['Status'] = dfs['Status']

    cols = ['FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched']
    df[cols] = dfs_only[cols]

    # Only plot screened genes
    df = df.loc[(df['Status'] == 'Screened'), :]

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

    # Only RNAi = True/FAlse
    if rnai is not None:
        df = df.loc[(df['RNAi']) == rnai, :]
    n_rows = len(df)

    codemap = {

    }
    code_label = {
        'A': 'Meiotic',
        'B': 'Post-meiotic',
        'C': 'Gametes',
        'D': 'Pre-meiotic',
        'E': 'General impairment of spermatogenesis',
        'F': 'Undetectable',
        'G': 'Unspecified ',
        'H': 'Non-germ cell autonomous'
    }
    code_color = {
        'A': '#d62728',
        'B': '#ce6dbd',
        'C': '#756bb1',
        'D': '#c7e9c0',
        'E': '#9edae5',
        'F': '#fdd0a2',
        'G': '#dadaeb',
        'H': '#bdbdbd'
    }

    #
    # Plot all data
    #
    print("Plotting")

    df = df.sort_values(['Status', 'mean fert-rate', 'gene'], ascending=[True, True, True]).reset_index()
    
    if rnai == 'Yes':
        figheight = 6
    elif rnai == 'No':
        figheight = 5
    else:
        figheight = 11
    fig = plt.figure(figsize=(2.2, figheight))
    # fig.suptitle('Core metazoan meiotic genes'.format(page, number_of_pages))

    gs = gridspec.GridSpec(nrows=1, ncols=12)
    ax_fert = plt.subplot(gs[:1, 0:11])
    ax_rnai = plt.subplot(gs[:1, 11])

    adjustable = 'datalim'
    aspect = 'auto'
    ax_fert.set(adjustable=adjustable, aspect=aspect, anchor='NE')
    ax_rnai.set(adjustable=adjustable, aspect=aspect, anchor='NE')

    rotation = 50
    s = 5
    marker = '_'
    n = len(df)
    yticks = list(np.arange(0, n, 50))
    #yticklabels = yticks[::-1]
    #
    # DM Expression Values
    #
    eb = ax_fert.errorbar(df['mean fert-rate'], range(0, len(df)), xerr=df['std fert-rate'], lw=0,
                          ecolor='#3182bd', elinewidth=0.6, capsize=0.6,
                          marker='.', markersize=4.0,
                          markeredgecolor='#3182bd', markeredgewidth=0.0,
                          markerfacecolor='#6baed6', markerfacecoloralt=None, zorder=5)
    ax_fert.axvline(0.75, color='#d62728', lw=1, zorder=6)
    ax_fert.set_xlabel('Fertility rate')
    ax_fert.set_ylabel('Core genes rank')
    ax_fert.set_xticks(np.linspace(0, 1, 5))
    # ax_fert.set_xticklabels([], fontsize='small', rotation=0)
    #ax_fert.set_yticks(yticks)
    #ax_fert.set_yticklabels(yticklabels)
    ax_fert.set_xlim(-0.04, 1.04)
    ax_fert.set_ylim(-1, len(df))
    ax_fert.grid(linewidth=0.5)
    ax_fert.set_xticklabels([], ha='right')
    ax_fert.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax_fert.tick_params(axis='both', labelsize='small')
    ax_fert.invert_yaxis()

    #
    # RNAi
    #
    data_rnai = df.loc[df['RNAi'] == 'Yes', 'RNAi']
    y = data_rnai.index
    x = np.zeros(len(y))
    c = '#17becf'

    sc_rnai = ax_rnai.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
    ax_rnai.set_xticks([0])
    ax_rnai.set_xticklabels([])
    ax_rnai.set_yticks(yticks)
    ax_rnai.tick_params(axis='y', which='major', length=1.5)
    ax_rnai.set_yticklabels([])
    ax_rnai.set_xlim(-0.2, 0.2)
    ax_rnai.set_ylim(-1, len(df))
    ax_rnai.set_xlabel('RNAi', rotation=90)
    # ax_rnai.grid(axis='y', linewidth=0.5)
    ax_rnai.invert_yaxis()

    if rnai == 'Yes':
        plt.subplots_adjust(left=0.25, right=0.96, bottom=0.09, top=0.98, wspace=0.9, hspace=0)
    elif rnai == 'No':
        plt.subplots_adjust(left=0.25, right=0.96, bottom=0.11, top=0.98, wspace=0.9, hspace=0)
    else:
        plt.subplots_adjust(left=0.25, right=0.96, bottom=0.050, top=0.99, wspace=0.9, hspace=0)
    #
    fig.savefig('images/img-core_DM_screened-vertical-rnai-{rnai:s}.pdf'.format(rnai=str(rnai)))
