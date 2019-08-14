# coding=utf-8
# Author: Rion B Correia
# Date: Jul 24, 2019
#
# Description: Plots HS and MM genes that are differentially expressed, calculated by 1_calc_diff_gene_exp.R
#
# Instructions:
#   For each species, two different comparisons are necessary:
#   Spermatocytes vs Spermatogonia (interested in genes upregulated in spermatocytes)
#   Spermatids vs Spermatocytes (interested in genes downregulated in spermatids)

#   Genes are considered upregulated if: log2FC >1 +  Pval < 0.01 + FDR≤0.05
#   Genes are considered downregulated if: log2FC <-1 +  Pval < 0.01 + FDR≤0.05
#
import math
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt


def plot_MA(df, file, title="plotMA", c_up='black', c_not='black', c_down='black'):
    s = 5
    lw = 0
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    maxFDR = 0.05
    minLogFC = math.log2(2)

    # Mask Data
    mask = ((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC))
    maskup = mask & (df['logFC'] >= 0)
    maskdown = mask & (df['logFC'] <= 0)
    masknot = ~mask
    # Filter Data
    dfUp = df.loc[maskup, :]
    dfNot = df.loc[masknot, :]
    dfDown = df.loc[maskdown, :]
    # Counts
    n_up, n_not, n_down = len(dfUp), len(dfNot), len(dfDown)
    # Plot
    ax.scatter(dfUp['logCPM'], dfUp['logFC'], c=c_up, s=s, zorder=2, lw=lw, alpha=.5)
    ax.scatter(dfNot['logCPM'], dfNot['logFC'], c=c_not, s=s / 3, zorder=1, lw=lw, alpha=.5)
    ax.scatter(dfDown['logCPM'], dfDown['logFC'], c=c_down, s=s, zorder=2, lw=lw, alpha=.5)
    # Draw a line at y=(-1,0,1)
    ax.axhline(y=-1, color='b', lw=1, linestyle='--')
    ax.axhline(y=0, color='gray', lw=1, linestyle='--')
    ax.axhline(y=+1, color='b', lw=1, linestyle='--')
    # Number of Selected Genes
    ax.text(x=0.96, y=0.93, s='{:,d}'.format(n_up), color=c_up, ha='right', va='center', transform=ax.transAxes, fontsize='large')
    ax.text(x=0.04, y=0.93, s='{:,d}'.format(n_not), color=c_not, ha='left', va='center', transform=ax.transAxes, fontsize='large')
    ax.text(x=0.96, y=0.07, s='{:,d}'.format(n_down), color=c_down, ha='right', va='center', transform=ax.transAxes, fontsize='large')
    # Labels
    ax.set_title(title)
    ax.set_ylabel('logFC')
    ax.set_xlabel('Average logCPM')
    # ax.grid()
    # Save
    plt.tight_layout()
    fig.savefig(file)


if __name__ == '__main__':

    #
    # [H]omo [S]apiens
    #
    # spermatocytes vs spermatogonia (interested in genes upregulated in spermatocytes)
    #
    df = pd.read_csv('results/HS-DGE_Cyte_vs_Gonia.csv', index_col=0, nrows=None)
    plot_MA(df=df, file='images/HS-DGE_UpCyte_vs_Gonia.pdf', title="HS (Up)Cyte vs Gonia", c_up='red', c_down='gray')

    #
    # Spermatids vs Spermatocytes (interested in genes downregulated in spermatids)
    #
    df = pd.read_csv('results/HS-DGE_Cyte_vs_Tid.csv', index_col=0, nrows=None)
    plot_MA(df=df, file='images/HS-DGE_DownCyte_vs_Tid.pdf', title="HS (Down)Gonia vs Tid", c_down='red', c_up='gray')

    #
    # MM
    #
    # spermatocytes vs spermatogonia (interested in genes upregulated in spermatocytes)
    #
    df = pd.read_csv('results/MM-DGE_Cyte_vs_Gonia.csv', index_col=0, nrows=None)
    plot_MA(df=df, file='images/MM-DGE_UpCyte_vs_Gonia.pdf', title="MM (Up)Cyte vs Gonia", c_up='red', c_down='gray')

    #
    # Spermatids vs Spermatocytes (interested in genes downregulated in spermatids)
    #
    df = pd.read_csv('results/MM-DGE_Cyte_vs_Tid.csv', index_col=0, nrows=None)
    plot_MA(df=df, file='images/MM-DGE_DownCyte_vs_Tid.pdf', title="MM (Down)Cyte vs Tid", c_down='red', c_up='gray')
