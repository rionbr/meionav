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


def plot_MA(df, file, title="plotMA", c_up5='black', c_up1='black', c_not='black', c_down5='black', c_down1='black'):
    s = 5
    lw = 0
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    minLogFC = math.log2(2)

    # Mask Data
    masklogfc = (df['logFC'].abs() >= minLogFC)
    mask5 = masklogfc & (df['FDR'] > 0.01) & (df['FDR'] <= 0.05)
    mask1 = masklogfc & (df['FDR'] <= 0.01) & (df['logFC'].abs() >= minLogFC)
    #
    maskup5 = mask5 & (df['logFC'] >= 0)
    maskup1 = mask1 & (df['logFC'] >= 0)
    maskdown5 = mask5 & (df['logFC'] <= 0)
    maskdown1 = mask1 & (df['logFC'] <= 0)
    masknot = ~(mask5 | mask1)
    # Filter Data
    dfUp5 = df.loc[maskup5, :]
    dfUp1 = df.loc[maskup1, :]
    dfNot = df.loc[masknot, :]
    dfDown5 = df.loc[maskdown5, :]
    dfDown1 = df.loc[maskdown1, :]
    # Counts
    n_up5, n_up1, n_not, n_down5, n_down1 = len(dfUp5), len(dfUp1), len(dfNot), len(dfDown5), len(dfDown1)
    # Plot
    ax.scatter(dfUp5['logCPM'], dfUp5['logFC'], c=c_up5, s=s, zorder=2, lw=lw, alpha=.5)
    ax.scatter(dfUp1['logCPM'], dfUp1['logFC'], c=c_up1, s=s, zorder=2, lw=lw, alpha=.5)
    ax.scatter(dfNot['logCPM'], dfNot['logFC'], c=c_not, s=s / 3, zorder=1, lw=lw, alpha=.5)
    ax.scatter(dfDown5['logCPM'], dfDown5['logFC'], c=c_down5, s=s, zorder=2, lw=lw, alpha=.5)
    ax.scatter(dfDown1['logCPM'], dfDown1['logFC'], c=c_down1, s=s, zorder=2, lw=lw, alpha=.5)
    # Draw a line at y=(-1,0,1)
    ax.axhline(y=-1, color='b', lw=1, linestyle='--')
    ax.axhline(y=0, color='gray', lw=1, linestyle='--')
    ax.axhline(y=+1, color='b', lw=1, linestyle='--')
    # Number of Selected Genes
    ax.text(x=0.96, y=0.93, s='{:,d}+{:,d}={:,d}'.format(n_up5, n_up1, n_up5 + n_up1), color=c_up5, ha='right', va='center', transform=ax.transAxes, fontsize='large')
    ax.text(x=0.04, y=0.93, s='{:,d}'.format(n_not), color=c_not, ha='left', va='center', transform=ax.transAxes, fontsize='large')
    ax.text(x=0.96, y=0.07, s='{:,d}+{:,d}={:,d}'.format(n_down5, n_down1, n_down5 + n_down1), color=c_down5, ha='right', va='center', transform=ax.transAxes, fontsize='large')
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
    plot_MA(df=df, file='images/HS-DGE_UpCyte_vs_Gonia.pdf', title="HS (Up)Cyte vs Gonia", c_up1='#d62728', c_up5='#ff9896', c_down1='gray', c_down5='gray')
    
    #
    # Spermatids vs Spermatocytes (interested in genes downregulated in spermatids)
    #
    df = pd.read_csv('results/HS-DGE_Cyte_vs_Tid.csv', index_col=0, nrows=None)
    plot_MA(df=df, file='images/HS-DGE_DownCyte_vs_Tid.pdf', title="HS (Down)Gonia vs Tid", c_down1='#1f77b4', c_down5='#aec7e8', c_up1='gray', c_up5='gray')

    #
    # MM
    #
    # spermatocytes vs spermatogonia (interested in genes upregulated in spermatocytes)
    #
    df = pd.read_csv('results/MM-DGE_Cyte_vs_Gonia.csv', index_col=0, nrows=None)
    plot_MA(df=df, file='images/MM-DGE_UpCyte_vs_Gonia.pdf', title="MM (Up)Cyte vs Gonia", c_up1='#d62728', c_up5='#ff9896', c_down1='gray', c_down5='gray')

    #
    # Spermatids vs Spermatocytes (interested in genes downregulated in spermatids)
    #
    df = pd.read_csv('results/MM-DGE_Cyte_vs_Tid.csv', index_col=0, nrows=None)
    plot_MA(df=df, file='images/MM-DGE_DownCyte_vs_Tid.pdf', title="MM (Down)Cyte vs Tid", c_down1='#1f77b4', c_down5='#aec7e8', c_up1='gray', c_up5='gray')
