# coding=utf-8
# Author: Rion B Correia
# Date: Jul 24, 2019
#
# Description: Plots HS and MM genes that are differentially expressed, calculated by 1_calc_diff_gene_exp.R
#
# Instructions:
#   For each species, two different comparisons are necessary:
#   Cytes vs Gonia (interested in genes upregulated in Cytes)
#   Cytes vs Tids (interested in genes downregulated in Cytes)
#   Genes are considered upregulated if: log2FC >1 +  Pval < 0.01 + FDR≤0.05
#   Genes are considered downregulated if: log2FC <-1 +  Pval < 0.01 + FDR≤0.05
#
import math
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
from utils import ensurePathExists


def plot_MA(df, core=[], pool=[], file='image.pdf', title="plotMA",
            c_up='#d62728', c_not='black', c_down='#1f77b4'
            ):
    s = 5
    lw = 0
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    minLogFC = math.log2(2)
    maxFDR = 0.05

    # Divide data into DGE Blocks
    dfU = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] >= 0), :]
    dfD = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] <= 0), :]
    dfN = df.loc[~df.index.isin(dfU.index.tolist() + dfD.index.tolist()), :]

    # Counts
    n_up, n_down, n_not = len(dfU), len(dfD), len(dfN)

    print("Up  : {rest:d} rest".format(rest=n_up))
    print("Down: {rest:d} rest".format(rest=n_down))
    print("Not : {rest:d} rest".format(rest=n_not))
    # Plot
    ax.scatter(dfU['logCPM'], dfU['logFC'], c=c_up, s=s, lw=lw, marker='o', zorder=3)
    ax.scatter(dfD['logCPM'], dfD['logFC'], c=c_down, s=s, lw=lw, marker='o', zorder=3)
    ax.scatter(dfN['logCPM'], dfN['logFC'], c=c_not, s=s / 3, lw=lw, marker='o', zorder=2)
    
    # Draw a line at y=(-1,0,1)
    ax.axhline(y=-1, color='b', lw=1, linestyle='--', zorder=5)
    ax.axhline(y=0, color='gray', lw=1, linestyle='--', zorder=5)
    ax.axhline(y=+1, color='b', lw=1, linestyle='--', zorder=5)

    ax.set_xlim(-1,18)
    ax.set_ylim(-15,15)


    # Labels
    ax.set_title(title)
    ax.set_ylabel('logFC')
    ax.set_xlabel('Average logCPM')
    # Layout
    plt.tight_layout()
    # Save
    
    ensurePathExists(file)
    fig.savefig(file)


if __name__ == '__main__':

    #
    # [H]omo [S]apiens
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])
    plot_MA(df=df, file='images/simpler/HS-DGE_Cyte_vs_Gonia.pdf', title="HS (Up)Cyte vs Gonia")

    #
    # Tid vs Cyte (interested in genes downregulated in Tid)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Tid_vs_Cyte.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])    
    plot_MA(df=df, file='images/simpler/HS-DGE_Tid_vs_Cyte.pdf', title="HS (Down)Tid vs Cyte")

    #
    # MM
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0)
    plot_MA(df=df, file='images/simpler/MM-DGE_Cyte_vs_Gonia.pdf', title="MM (Up)Cyte vs Gonia")
    #
    # Cytes vs Tids (interested in genes downregulated in Tid)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Tid_vs_Cyte.csv', index_col=0)
    plot_MA(df=df, file='images/simpler/MM-DGE_Tid_vs_Cyte.pdf', title="MM (Down)Tid vs Cyte")

    #
    # DM
    #
    # Middle vs Apical (interested in genes upregulated in Middle)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Apical.csv', index_col=0)
    plot_MA(df=df, file='images/simpler/DM-DGE_Middle_vs_Apical.pdf', title="DM (Up)Middle vs Apical")

    #
    # Basal vs Middle (interested in genes downregulated in Basal)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Basal_vs_Middle.csv', index_col=0)
    plot_MA(df=df, file='images/simpler/DM-DGE_Basal_vs_Middle.pdf', title="DM (Down)Basal vs Middle")
