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


def plot_MA(df, core=[], pool=[], file='image.pdf', title="plotMA",
            x_up=0, y_up=0,
            x_not=0, y_not=0,
            x_down=0, y_down=0,
            c_core='black', c_pool='black', c_up='black', c_not='black', c_down='black',
            m_core='P', m_pool='X', m_up='o', m_not='o', m_down='o'
            ):
    s = 5
    lw = 0
    alpha = 0.8
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    minLogFC = math.log2(2)
    maxFDR = 0.05

    # Divide data into DGE Blocks
    dfU = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] >= 0), :]
    dfD = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] <= 0), :]
    dfN = df.loc[~df.index.isin(dfU.index.tolist() + dfD.index.tolist()), :]

    # Filter Data
    dfUC = dfU.loc[(dfU.index.isin(core)), :]
    dfUP = dfU.loc[(dfU.index.isin(pool)) & ~(dfU.index.isin(core)), :]
    dfUN = dfU.loc[~dfU.index.isin(dfUC.index.tolist() + dfUP.index.tolist()), :]

    dfDC = dfD.loc[(dfD.index.isin(core)), :]
    dfDP = dfD.loc[(dfD.index.isin(pool)) & ~(dfD.index.isin(core)), :]
    dfDN = dfD.loc[~dfD.index.isin(dfDC.index.tolist() + dfDP.index.tolist()), :]

    dfNC = dfN.loc[(dfN.index.isin(core)), :]
    dfNP = dfN.loc[(dfN.index.isin(pool)) & ~(dfN.index.isin(core)), :]
    dfNN = dfN.loc[~dfN.index.isin(dfNC.index.tolist() + dfNP.index.tolist()), :]

    # Sanity Check
    if not set(df.index.tolist()) == set(
            dfUC.index.tolist() + dfUP.index.tolist() + dfUN.index.tolist() +
            dfDC.index.tolist() + dfDP.index.tolist() + dfDN.index.tolist() +
            dfNC.index.tolist() + dfNP.index.tolist() + dfNN.index.tolist()
        ):
        raise ValueError('Sanity check has failed!')
    # Counts
    n_up_core, n_up_pool, n_up_rest = len(dfUC), len(dfUP), len(dfUN)
    n_down_core, n_down_pool, n_down_rest = len(dfDC), len(dfDP), len(dfDN)
    n_not_core, n_not_pool, n_not_rest = len(dfNC), len(dfNP), len(dfNN)

    print("Up  : {core:d} core, {pool:d} pool, {rest:d} rest".format(core=n_up_core, pool=n_up_pool, rest=n_up_rest))
    print("Down: {core:d} core, {pool:d} pool, {rest:d} rest".format(core=n_down_core, pool=n_down_pool, rest=n_down_rest))
    print("Not : {core:d} core, {pool:d} pool, {rest:d} rest".format(core=n_not_core, pool=n_not_pool, rest=n_not_rest))
    # Plot
    ax.scatter(dfUC['logCPM'], dfUC['logFC'], c=c_core, s=s * 2, lw=lw, alpha=alpha, marker=m_core, zorder=4)
    ax.scatter(dfUP['logCPM'], dfUP['logFC'], c=c_pool, s=s * 2, lw=lw, alpha=alpha, marker=m_pool, zorder=3)
    ax.scatter(dfUN['logCPM'], dfUN['logFC'], c=c_up, s=s, lw=lw, alpha=alpha, marker=m_up, zorder=2)
    #
    ax.scatter(dfDC['logCPM'], dfDC['logFC'], c=c_core, s=s * 2, lw=lw, alpha=alpha, marker=m_core, zorder=4)
    ax.scatter(dfDP['logCPM'], dfDP['logFC'], c=c_pool, s=s * 2, lw=lw, alpha=alpha, marker=m_pool, zorder=3)
    ax.scatter(dfDN['logCPM'], dfDN['logFC'], c=c_down, s=s, lw=lw, alpha=alpha, marker=m_down, zorder=2)

    ax.scatter(dfNC['logCPM'], dfNC['logFC'], c=c_core, s=s * 2, lw=lw, alpha=alpha, marker=m_core, zorder=4)
    ax.scatter(dfNP['logCPM'], dfNP['logFC'], c=c_pool, s=s * 2, lw=lw, alpha=alpha, marker=m_pool, zorder=3)
    ax.scatter(dfNN['logCPM'], dfNN['logFC'], c=c_not, s=s / 3, lw=lw, alpha=alpha, marker=m_not, zorder=2)

    # Draw a line at y=(-1,0,1)
    ax.axhline(y=-1, color='b', lw=1, linestyle='--')
    ax.axhline(y=0, color='gray', lw=1, linestyle='--')
    ax.axhline(y=+1, color='b', lw=1, linestyle='--')

    ax.set_xlim(-1, 18)
    ax.set_ylim(-15, 15)

    # Number of Selected Genes
    """
    # Up
    strings, colors = [], []
    n_up_total = 0
    if n_up_core > 0:
        strings.extend(['{:,d}'.format(n_up_core), '+'])
        colors.extend([c_core, 'black'])
        n_up_total += n_up_core
    if n_up_pool > 0:
        strings.extend(['{:,d}'.format(n_up_pool), '+'])
        colors.extend([c_pool, 'black'])
        n_up_total += n_up_pool
    n_up_total += n_up_rest
    strings.extend(['{:,d}'.format(n_up_rest), '=', '{:,d}'.format(n_up_total)])
    colors.extend([c_up, 'black', 'black'])
    rainbow_text(ax=ax, x=8, y=13, strings=strings, colors=colors, ha='left', va='center', fontsize='large')

    # Down
    strings, colors = [], []
    n_down_total = 0
    if n_down_core > 0:
        strings.extend(['{:,d}'.format(n_down_core), '+'])
        colors.extend([c_core, 'black'])
        n_down_total += n_down_core
    if n_down_pool > 0:
        strings.extend(['{:,d}'.format(n_down_pool), '+'])
        colors.extend([c_pool, 'black'])
        n_down_total += n_down_pool
    n_down_total += n_down_rest
    strings.extend(['{:,d}'.format(n_down_rest), '=', '{:,d}'.format(n_down_total)])
    colors.extend([c_down, 'black', 'black'])
    rainbow_text(ax=ax, x=8, y=-13, strings=strings, colors=colors, ha='left', va='center', fontsize='large')

    # Not
    strings, colors = [], []
    n_not_total = 0
    if n_not_core > 0:
        strings.extend(['{:,d}'.format(n_not_core), '+'])
        colors.extend([c_core, 'black'])
        n_not_total += n_not_core
    if n_not_pool > 0:
        strings.extend(['{:,d}'.format(n_not_pool), '+'])
        colors.extend([c_pool, 'black'])
        n_not_total += n_not_pool
    n_not_total += n_not_rest
    strings.extend(['{:,d}'.format(n_not_rest), '=', '{:,d}'.format(n_not_total)])
    colors.extend([c_not, 'black', 'black'])
    rainbow_text(ax=ax, x=10, y=-3, strings=strings, colors=colors, ha='left', va='center', fontsize='large')
    """
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

    core_pipeline = 'all3-conserved'
    pool_pipeline = 'all3-pooling-DM'
    #
    # [H]omo [S]apiens
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0)
    core = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] >= 1) & (dfC['FDR_CyteGonia'] <= 0.05), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0)
    pool = dfP.loc[(dfP['Cyte_vs_Gonia'] == True) & (dfP['logFC_CyteGonia'] > 1) & (dfP['FDR_CyteGonia'] < 0.05), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/HS-DGE_UpCyte_vs_Gonia.pdf', title="HS (Up)Cyte vs Gonia",
            x_up=10, y_up=5, x_not=15, y_not=0, x_down=15, y_down=-7,
            c_core='#d62728', c_pool='#9467bd', c_up='#ff9896', c_down='gray'
            )
    #
    # Tid vs Cyte (interested in genes downregulated in Tid)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Tid_vs_Cyte.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0)
    core = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] <= -1) & (dfC['FDR_TidCyte'] <= 0.05), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0)
    pool = dfP.loc[(dfP['Tid_vs_Cyte'] == True) & (dfP['logFC_TidCyte'] <= -1) & (dfP['FDR_TidCyte'] <= 0.05), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/HS-DGE_DownTid_vs_Cyte.pdf', title="HS (Down)Tid vs Cyte",
            x_up=14, y_up=6, x_not=14, y_not=0, x_down=9, y_down=-13,
            c_core='#1f77b4', c_pool='#9467bd', c_up='gray', c_down='#aec7e8'
            )
    #
    # MM
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0)
    dfC = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0)
    core = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] >= 1) & (dfC['FDR_CyteGonia'] <= 0.05), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0)
    pool = dfP.loc[(dfP['Cyte_vs_Gonia'] == True) & (dfP['logFC_CyteGonia'] >= 1) & (dfP['FDR_CyteGonia'] <= 0.05), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/MM-DGE_UpCyte_vs_Gonia.pdf', title="MM (Up)Cyte vs Gonia",
            x_up=1, y_up=9, x_not=11, y_not=0, x_down=1, y_down=-9,
            c_core='#d62728', c_pool='#9467bd', c_up='#ff9896', c_down='gray'
            )
    #
    # Cytes vs Tids (interested in genes downregulated in Tid)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Tid_vs_Cyte.csv', index_col=0)
    dfC = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0)
    core = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] <= -1) & (dfC['FDR_TidCyte'] <= 0.05), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0)
    pool = dfP.loc[(dfP['Tid_vs_Cyte'] == True) & (dfP['logFC_TidCyte'] <= -1) & (dfP['FDR_TidCyte'] <= 0.05), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/MM-DGE_DownTid_vs_Cyte.pdf', title="MM (Down)Tid vs Cyte",
            x_up=11, y_up=4, x_not=11, y_not=0, x_down=7, y_down=-6,
            c_core='#1f77b4', c_pool='#9467bd', c_down='#aec7e8', c_up='gray'
            )
    #
    # DM
    #
    # Middle vs Apical (interested in genes upregulated in Middle)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Apical.csv', index_col=0)
    dfC = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0)
    core = dfC.loc[(dfC['Middle_vs_Apical'] == True) & (dfC['logFC_MiddleApical'] >= 1) & (dfC['FDR_MiddleApical'] <= 0.05), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0)
    pool = dfP.loc[(dfP['Middle_vs_Apical'] == True) & (dfP['logCPM_MiddleApical'] >= 1), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/DM-DGE_UpMiddle_vs_Apical.pdf', title="DM (Up)Middle vs Apical",
            x_up=8, y_up=4, x_not=12, y_not=-0.40, x_down=12, y_down=-5,
            c_core='#d62728', c_pool='#9467bd',
            c_up='#ff9896', c_down='gray'
            )
    #
    # Basal vs Middle (interested in genes downregulated in Basal)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Basal_vs_Middle.csv', index_col=0)
    dfC = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0)
    core = dfC.loc[(dfC['Basal_vs_Middle'] == True) & (dfC['logFC_BasalMiddle'] <= -1) & (dfC['FDR_BasalMiddle'] <= 0.05), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0)
    pool = dfP.loc[(dfP['Basal_vs_Middle'] == True) & (dfP['logCPM_BasalMiddle'] >= 1), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/DM-DGE_DownBasal_vs_Middle.pdf', title="DM (Down)Basal vs Middle",
            x_up=12.5, y_up=4, x_not=12, y_not=0.375, x_down=9, y_down=-6,
            c_core='#1f77b4', c_pool='#9467bd',
            c_down='#aec7e8', c_up='gray'
            )
