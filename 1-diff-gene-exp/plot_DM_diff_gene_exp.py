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
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt


def plot_MA(df, file, title="plotMA", c_p='black', c_a='black'):
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))
    df.sort_values('Ave signal', ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    # Mask Data
    maskA = (df['Combined detection signal'] == 'A')
    maskP = (df['Combined detection signal'] == 'P')
    # Filter Data
    dfA = df.loc[maskA, :]
    dfP = df.loc[maskP, :]
    # Counts
    n_a, n_p = len(dfA), len(dfP)
    # Plot
    # ax.scatter(np.random.rand(n_genes), np.log(df['Ave signal']), c=c_up, s=s, zorder=2)
    ax.scatter(dfA.index, dfA['Ave signal'], c=c_a, s=s, zorder=2)
    ax.scatter(dfP.index, dfP['Ave signal'], c=c_p, s=s, zorder=2)
    # Draw a line at y=(-1,0,1)
    # ax.axhline(y=-1, color='b', linestyle='--')
    # ax.axhline(y=0, color='yellow', linestyle='--')
    # ax.axhline(y=+1, color='b', linestyle='--')
    # Number of Selected Genes
    ax.text(x=0.96, y=0.87, s='A={:,d}'.format(n_a), color=c_a, ha='right', va='center', transform=ax.transAxes, fontsize='large')
    ax.text(x=0.96, y=0.93, s='P={:,d}'.format(n_p), color=c_p, ha='right', va='center', transform=ax.transAxes, fontsize='large')
    # Labels
    ax.set_title(title)
    ax.set_ylabel('Average signal')
    ax.set_xlabel('Gene rank')
    #
    ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(0.1,10e5)
    ax.grid()
    # Save
    plt.tight_layout()
    fig.savefig(file)


if __name__ == '__main__':

    #
    # Init Plot Params
    #
    s = 6

    #
    # [D]rosophila [M]elanogaster
    #
    # Apical Testis
    #
    df = pd.read_csv('data/DM_ApicalTestis_Exp_Pcall.csv', index_col=0, na_values='---', nrows=None)
    df = df.dropna(subset=['Ensembl'])
    df = df.drop_duplicates(subset=['Ensembl'])
    plot_MA(df=df, file='images/DM-apical.pdf', title="DM apical testis", c_p='red', c_a='blue')

    #
    # Mid Testis
    #
    df = pd.read_csv('data/DM_MidTestis_Exp_Pcall.csv', index_col=0, na_values='---', nrows=None)
    df = df.dropna(subset=['Ensembl'])
    df = df.drop_duplicates(subset=['Ensembl'])
    plot_MA(df=df, file='images/DM-mid.pdf', title="DM mid testis", c_p='red', c_a='blue')
