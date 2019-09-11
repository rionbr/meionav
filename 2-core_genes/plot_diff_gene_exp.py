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
            x_up=0, y_up=0,
            x_not=0, y_not=0,
            x_down=0, y_down=0,
            c_up='black', c_up_core='black',
            c_not='black',
            c_down='black', c_down_core='black'):
    s = 5
    lw = 0
    alpha = 0.8
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    minLogFC = math.log2(2)

    # Mask Data
    mask_core = (df.index.isin(core))
    mask_logfc = (df['logFC'].abs() >= minLogFC)
    mask_fdr = (df['FDR'] <= 0.05)
    #
    mask_up = ~mask_core & mask_logfc & mask_fdr & (df['logFC'] >= 0)
    mask_up_core = mask_core & mask_logfc & mask_fdr & (df['logFC'] >= 0)
    mask_down = ~mask_core & mask_logfc & mask_fdr & (df['logFC'] <= 0)
    mask_down_core = mask_core & mask_logfc & mask_fdr & (df['logFC'] <= 0)
    mask_not = ~(mask_up | mask_up_core | mask_down | mask_down_core)
    # Filter Data
    dfUp = df.loc[mask_up, :]
    dfUpC = df.loc[mask_up_core, :]
    dfNot = df.loc[mask_not, :]
    dfDown = df.loc[mask_down, :]
    dfDownC = df.loc[mask_down_core, :]
    # Counts
    n_up, n_up_core, n_not, n_down, n_down_core = len(dfUp), len(dfUpC), len(dfNot), len(dfDown), len(dfDownC)
    # Plot
    ax.scatter(dfUp['logCPM'], dfUp['logFC'], c=c_up, s=s, lw=lw, alpha=alpha, zorder=2)
    ax.scatter(dfUpC['logCPM'], dfUpC['logFC'], c=c_up_core, s=s * 2, lw=lw, alpha=alpha, zorder=3)
    ax.scatter(dfNot['logCPM'], dfNot['logFC'], c=c_not, s=s / 3, lw=lw, alpha=alpha, zorder=1)
    ax.scatter(dfDown['logCPM'], dfDown['logFC'], c=c_down, s=s, lw=lw, alpha=alpha, zorder=2)
    ax.scatter(dfDownC['logCPM'], dfDownC['logFC'], c=c_down_core, s=s * 2, lw=lw, alpha=alpha, zorder=3)
    # Draw a line at y=(-1,0,1)
    ax.axhline(y=-1, color='b', lw=1, linestyle='--')
    ax.axhline(y=0, color='gray', lw=1, linestyle='--')
    ax.axhline(y=+1, color='b', lw=1, linestyle='--')
    
    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_ylim())
    
    # Number of Selected Genes
    # Up
    if n_up_core > 0:
        strings = [
            '{:,d}'.format(n_up),
            '+',
            '{:,d}'.format(n_up_core),
            '=',
            '{:,d}'.format(n_up + n_up_core)]
        colors = [c_up, 'black', c_up_core, 'black', 'black']
    else:
        strings = ['{:,d}'.format(n_up)]
        colors = ['black']
    rainbow_text(ax=ax, x=x_up, y=y_up, strings=strings, colors=colors, ha='left', va='center', fontsize='large')
    # Not
    ax.text(x=x_not, y=y_not, s='{:,d}'.format(n_not), color=c_not, ha='left', va='center', fontsize='large')
    # Down
    if n_down_core > 0:
        strings = [
            '{:,d}'.format(n_down),
            '+',
            '{:,d}'.format(n_down_core),
            '=',
            '{:,d}'.format(n_down + n_down_core)]
        colors = [c_down, 'black', c_down_core, 'black', 'black']
    else:
        strings = ['{:,d}'.format(n_down)]
        colors = ['black']
    rainbow_text(ax=ax, x=x_down, y=y_down, strings=strings, colors=colors, ha='left', va='center', fontsize='large')
    
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

    pipeline = 'all3-conserved-FDRp05'
    #
    # [H]omo [S]apiens
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0, nrows=None)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['UpCyte_vs_Gonia'] == True) & (dfC['logFC_up'] > 0), :].index.tolist()
    
    plot_MA(df=df, core=core, file='images/{pipeline:s}/HS-DGE_UpCyte_vs_Gonia.pdf'.format(pipeline=pipeline), title="HS (Up)Cyte vs Gonia",
        x_up=10, y_up=5, x_not=15, y_not=0, x_down=15, y_down=-7,
        c_up_core='#d62728', c_up='#ff9896', c_down='gray'
        )
    #
    # Cytes vs Tids (interested in genes downregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Tid.csv', index_col=0, nrows=None)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['DownCyte_vs_Tid'] == True) & (dfC['logFC_down'] < 0), :].index.tolist()

    plot_MA(df=df, core=core, file='images/{pipeline:s}/HS-DGE_DownCyte_vs_Tid.pdf'.format(pipeline=pipeline), title="HS (Down)Gonia vs Tid",
        x_up=14, y_up=6, x_not=14, y_not=0, x_down=9, y_down=-13,
        c_down_core='#1f77b4', c_down='#aec7e8', c_up='gray'
        )

    #
    # MM
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['UpCyte_vs_Gonia'] == True) & (dfC['logFC_up'] > 0), :].index.tolist()
    
    plot_MA(df=df, core=core, file='images/{pipeline:s}/MM-DGE_UpCyte_vs_Gonia.pdf'.format(pipeline=pipeline), title="MM (Up)Cyte vs Gonia",
        x_up=1, y_up=9, x_not=11, y_not=0, x_down=1, y_down=-9,
        c_up_core='#d62728', c_up='#ff9896', c_down='gray'
        )
    #
    # Cytes vs Tids (interested in genes downregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Tid.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['DownCyte_vs_Tid'] == True) & (dfC['logFC_down'] < 0), :].index.tolist()

    plot_MA(df=df, core=core, file='images/{pipeline:s}/MM-DGE_DownCyte_vs_Tid.pdf'.format(pipeline=pipeline), title="MM (Down)Cyte vs Tid",
        x_up=11, y_up=4, x_not=11, y_not=0, x_down=7, y_down=-6,
        c_down_core='#1f77b4', c_down='#aec7e8', c_up='gray'
        )

    #
    # DM
    #
    # Middle vs Apical (interested in genes upregulated in Middle)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Apical.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['UpMiddle_vs_Apical'] == True) & (dfC['logFC_up'] > 0), :].index.tolist()
    
    plot_MA(df=df, core=core, file='images/{pipeline:s}/DM-DGE_UpMiddle_vs_Apical.pdf'.format(pipeline=pipeline), title="DM (Up)Middle vs Apical",
        x_up=8, y_up=4, x_not=12, y_not=-0.40, x_down=12, y_down=-5,
        c_up_core='#d62728', c_up='#ff9896', c_down='gray'
        )
    #
    # Middle vs Basal (interested in genes downregulated in Middle)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Basal.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['DownMiddle_vs_Basal'] == True) & (dfC['logFC_down'] < 0), :].index.tolist()

    plot_MA(df=df, core=core, file='images/{pipeline:s}/DM-DGE_DownMiddle_vs_Basal.pdf'.format(pipeline=pipeline), title="DM (Down)Middle vs Basal",
        x_up=12.5, y_up=4, x_not=12, y_not=0.375, x_down=9, y_down=-6,
        c_down_core='#1f77b4', c_down='#aec7e8', c_up='gray'
        )