# coding=utf-8
# Author: Rion B Correia
# Date: Jul 24, 2019
#
# Description: Plots DM genes that are differentially expressed.
#
# Instructions:
#   Raw data for DM: apical and medial testis gene expression profiles.
#   - y-axis: average signal
#   - x-axis: genes (Ensembl), colour-coded by either a Present or Absent call (combined detection signal).
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
#mpl.use('agg')
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from matplotlib import transforms


def rainbow_text(x, y, strings, colors, orientation='horizontal',
                 ax=None, **kwargs):
    t = ax.transData
    canvas = ax.figure.canvas

    for s, c in zip(strings, colors):
        text = ax.text(x, y, s, color=c, **kwargs)

        # Need to draw to update the text position.
        text.draw(canvas.get_renderer())
        bbox = text.get_window_extent() # in display units
        t = transforms.offset_copy(text.get_transform(), x=bbox.width, units='dots')
        bbox = bbox.inverse_transformed(ax.transData) # in data units
        x += bbox.width

def plot_MA(df, core=[], file='image.pdf', title="plotMA",
            x_p=0, y_p=0, x_a=0, y_a=0,
            c_p_core='black', c_p='black', c_a='black'):
    s = 8
    lw = 0
    alpha = 0.8
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
    df = df.sort_values('avg-signal', ascending=False).reset_index()
    # Mask Data
    mask_p = (df['A/P'] == 'P')
    mask_core = (df['Gene stable ID'].isin(core))
    mask_a = (df['A/P'] == 'A')
    # Filter Data
    df_p = df.loc[mask_p & ~mask_core, :]
    df_p_core = df.loc[mask_p & mask_core, :]
    df_a = df.loc[mask_a, :]
    # Counts
    n_p, n_p_core, n_a = len(df_p), len(df_p_core), len(df_a)
    # Plot
    ax.scatter(df_p_core.index, df_p_core['avg-signal'], c=c_p_core, s=s, lw=lw, alpha=alpha, zorder=3)
    ax.scatter(df_p.index, df_p['avg-signal'], c=c_p, s=s*2, lw=lw, alpha=alpha, zorder=2)
    ax.scatter(df_a.index, df_a['avg-signal'], c=c_a, s=s*4, lw=lw, alpha=alpha, zorder=1)

    # Number of Selected Genes
    strings = ['P=', '{:,d}'.format(n_p_core), '+', '{:,d}'.format(n_p), '=', '{:,d}'.format(n_p + n_p_core)]
    colors = ['black', c_p_core, 'black', c_p, 'black', 'black']
    rainbow_text(x=x_p, y=y_p, strings=strings, colors=colors, ax=ax, fontsize='large')
    ax.text(x=x_a, y=y_a, s='A={:,d}'.format(n_a), color=c_a, ha='left', va='center', fontsize='large')
    
    # Labels
    ax.set_title(title)
    ax.set_ylabel('Average signal')
    ax.set_xlabel('Gene rank')
    #
    ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(0.1,10e5)
    ylimmin,ylimmax = ax.get_ylim()
    ax.set_ylim(.3,ylimmax)
    
    ax.grid()
    # Save
    plt.tight_layout()
    fig.savefig(file)


if __name__ == '__main__':

    #
    # [D]rosophila [M]elanogaster
    dfC = pd.read_csv('results/DM/core_DM_meiotic_genes.csv', index_col=0, nrows=None)
    #
    # Apical Testis
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM_DGE_ApicalTestis.csv', nrows=None)
    df = df.drop_duplicates(subset=['Gene stable ID']).set_index('Gene stable ID') 
    core = dfC.loc[(dfC['DM_up'] == True), :].index.tolist()

    plot_MA(df=df, core=core, file='images/DM-apical.pdf', title="DM apical testis",
            x_p=2900, y_p=3000, x_a=10500, y_a=200,
            c_p_core='#d62728', c_p='#ff9896', c_a='gray')

    #
    # Mid Testis
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM_DGE_MidTestis.csv', nrows=None)
    df = df.drop_duplicates(subset=['Gene stable ID']).set_index('Gene stable ID')
    core = dfC.loc[(dfC['DM_down'] == True), :].index.tolist()
    plot_MA(df=df, core=core, file='images/DM-mid.pdf', title="DM mid testis",
            x_p=2900, y_p=3000, x_a=10500, y_a=200,
            c_p_core='#1f77b4', c_p='#aec7e8', c_a='gray')
