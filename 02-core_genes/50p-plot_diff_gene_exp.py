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
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
from utils import ensurePathExists


def rainbow_text(x, y, strings, colors, orientation='horizontal',
                 ax=None, **kwargs):
    canvas = ax.figure.canvas

    for s, c in zip(strings, colors):
        text = ax.text(x, y, s, color=c, **kwargs)

        # Need to draw to update the text position.
        text.draw(canvas.get_renderer())
        bbox = text.get_window_extent()  # in display units
        # transform to data units
        bbox = mpl.transforms.Bbox(ax.transData.inverted().transform(bbox))
        x += bbox.width


def plot_MA(df, highlight=[], pool=[], file='image.pdf', title="plotMA",
            annotate_plot=False,
            c_highlight='gray', c_up='gray', c_not='black', c_down='gray',
            m_highlight='P', m_pool='X', m_up='o', m_not='o', m_down='o'
            ):
    s = 5
    lw = 0
    alpha = 0.8
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))

    minLogFC = math.log2(2)
    maxFDR = 0.05

    # Divide data into DGE Blocks
    dfU = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] >= 0), :]
    dfD = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] <= 0), :]
    dfN = df.loc[~df.index.isin(dfU.index.tolist() + dfD.index.tolist()), :]

    # Filter Data
    dfUH = dfU.loc[(dfU.index.isin(highlight)), :]
    dfUN = dfU.loc[~dfU.index.isin(dfUH.index.tolist()), :]

    dfDH = dfD.loc[(dfD.index.isin(highlight)), :]
    dfDN = dfD.loc[~dfD.index.isin(dfDH.index.tolist()), :]

    dfNH = dfN.loc[(dfN.index.isin(highlight)), :]
    dfNN = dfN.loc[~dfN.index.isin(dfNH.index.tolist()), :]

    # Sanity Check
    if not set(df.index.tolist()) == set(
            dfUH.index.tolist() + dfUN.index.tolist() +
            dfDH.index.tolist() + dfDN.index.tolist() +
            dfNH.index.tolist() + dfNN.index.tolist()):
        raise ValueError('Sanity check has failed!')
    # Counts
    n_up_core, n_up_rest = len(dfUH), len(dfUN)
    n_down_core, n_down_rest = len(dfDH), len(dfDN)
    n_not_core, n_not_rest = len(dfNH), len(dfNN)

    print("Up  : {core:d} core, {rest:d} rest".format(core=n_up_core, rest=n_up_rest))
    print("Down: {core:d} core, {rest:d} rest".format(core=n_down_core, rest=n_down_rest))
    print("Not : {core:d} core, {rest:d} rest".format(core=n_not_core, rest=n_not_rest))
    # Plot
    ax.scatter(dfUH['logCPM'], dfUH['logFC'], c=c_highlight, s=s * 2, lw=lw, alpha=alpha, marker=m_highlight, zorder=4, rasterized=True)
    ax.scatter(dfUN['logCPM'], dfUN['logFC'], c=c_up, s=s, lw=lw, alpha=alpha, marker=m_up, zorder=2, rasterized=True)
    #
    ax.scatter(dfDH['logCPM'], dfDH['logFC'], c=c_highlight, s=s * 2, lw=lw, alpha=alpha, marker=m_highlight, zorder=4, rasterized=True)
    ax.scatter(dfDN['logCPM'], dfDN['logFC'], c=c_down, s=s, lw=lw, alpha=alpha, marker=m_down, zorder=2, rasterized=True)

    ax.scatter(dfNH['logCPM'], dfNH['logFC'], c=c_highlight, s=s * 2, lw=lw, alpha=alpha, marker=m_highlight, zorder=4, rasterized=True)
    ax.scatter(dfNN['logCPM'], dfNN['logFC'], c=c_not, s=s / 3, lw=lw, alpha=alpha, marker=m_not, zorder=2, rasterized=True)

    # Draw a line at y=(-1,0,1)
    ax.axhline(y=-1, color='b', lw=1, linestyle='--')
    ax.axhline(y=0, color='gray', lw=1, linestyle='--')
    ax.axhline(y=+1, color='b', lw=1, linestyle='--')

    ax.set_xlim(-1, 18)
    ax.set_ylim(-15, 15)

    # Number of Selected Genes
    if annotate_plot:
        # Up
        strings, colors = [], []
        n_up_total = 0
        if n_up_core > 0:
            strings.extend(['{:,d}'.format(n_up_core), '+'])
            colors.extend([c_highlight, 'black'])
            n_up_total += n_up_core
        n_up_total += n_up_rest
        strings.extend(['{:,d}'.format(n_up_rest), '=', '{:,d}'.format(n_up_total)])
        colors.extend([c_up, 'black', 'black'])
        rainbow_text(ax=ax, x=9, y=13, strings=strings, colors=colors, ha='left', va='center', fontsize='large')

        # Not
        strings, colors = [], []
        n_not_total = 0
        if n_not_core > 0:
            strings.extend(['{:,d}'.format(n_not_core), '+'])
            colors.extend([c_highlight, 'black'])
            n_not_total += n_not_core
        n_not_total += n_not_rest
        strings.extend(['{:,d}'.format(n_not_rest), '=', '{:,d}'.format(n_not_total)])
        colors.extend([c_not, 'black', 'black'])
        rainbow_text(ax=ax, x=10, y=-3, strings=strings, colors=colors, ha='left', va='center', fontsize='large')

        # Down
        strings, colors = [], []
        n_down_total = 0
        if n_down_core > 0:
            strings.extend(['{:,d}'.format(n_down_core), '+'])
            colors.extend([c_highlight, 'black'])
            n_down_total += n_down_core
        n_down_total += n_down_rest
        strings.extend(['{:,d}'.format(n_down_rest), '=', '{:,d}'.format(n_down_total)])
        colors.extend([c_down, 'black', 'black'])
        rainbow_text(ax=ax, x=9, y=-13, strings=strings, colors=colors, ha='left', va='center', fontsize='large')

    # Labels
    ax.set_title(title)
    ax.set_ylabel('logFC')
    ax.set_xlabel('Average logCPM')
    # Layout
    plt.subplots_adjust(left=0.17, bottom=0.17, right=0.97, top=0.90)
    #plt.tight_layout()
    # Save

    ensurePathExists(file)
    fig.savefig(file, dpi=300)


if __name__ == '__main__':

    #
    # Pipeline = Mammals
    #
    pipeline = 'mammals'
    c_highlight = '#9467bd'
    #
    # [H]omo [S]apiens
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../01-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    mammals = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] >= 1) & (dfC['FDR_CyteGonia'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/annotated/HS-DGE_UpCyte_vs_Gonia-ann.pdf', title="Human (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=True)
    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/clean/HS-DGE_UpCyte_vs_Gonia-cln.pdf', title="Human (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=False)

    #
    # Tid vs Cyte (interested in genes downregulated in Tid)
    #
    df = pd.read_csv('../01-diff-gene-exp/results/HS/HS-DGE_Tid_vs_Cyte.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    mammals = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] <= -1) & (dfC['FDR_TidCyte'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/annotated/HS-DGE_DownTid_vs_Cyte-ann.pdf', title="Human (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=True)
    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/clean/HS-DGE_DownTid_vs_Cyte-clean.pdf', title="Human (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=False)

    #
    # MM
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../01-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0)
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    mammals = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] >= 1) & (dfC['FDR_CyteGonia'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/annotated/MM-DGE_UpCyte_vs_Gonia-ann.pdf', title="Mouse (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=True)
    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/clean/MM-DGE_UpCyte_vs_Gonia-cln.pdf', title="Mouse (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=False)

    #
    # Cytes vs Tids (interested in genes downregulated in Tid)
    #
    df = pd.read_csv('../01-diff-gene-exp/results/MM/MM-DGE_Tid_vs_Cyte.csv', index_col=0)
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    mammals = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] <= -1) & (dfC['FDR_TidCyte'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/annotated/MM-DGE_DownTid_vs_Cyte-ann.pdf', title="Mouse (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_down='gray', c_up='#aec7e8', annotate_plot=True)
    plot_MA(df=df, highlight=mammals, file='images/pipeline-mammals/clean/MM-DGE_DownTid_vs_Cyte-cln.pdf', title="Mouse (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_down='gray', c_up='#aec7e8', annotate_plot=False)

    #
    # Pipeline = Core
    #
    pipeline = 'core'
    c_highlight = '#2ca02c'

    #
    # [H]omo [S]apiens
    #
    # Cytes vs Gonia
    #
    df = pd.read_csv('../01-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    core = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] >= 1) & (dfC['FDR_CyteGonia'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=core, file='images/pipeline-core/annotated/HS-DGE_UpCyte_vs_Gonia-ann.pdf', title="Human (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=True)
    plot_MA(df=df, highlight=core, file='images/pipeline-core/clean/HS-DGE_UpCyte_vs_Gonia-cln.pdf', title="Human (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=False)

    #
    # Tid vs Cyte
    #
    df = pd.read_csv('../01-diff-gene-exp/results/HS/HS-DGE_Tid_vs_Cyte.csv', index_col=0)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    core = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] <= -1) & (dfC['FDR_TidCyte'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=core, file='images/pipeline-core/annotated/HS-DGE_DownTid_vs_Cyte-ann.pdf', title="Human (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=True)
    plot_MA(df=df, highlight=core, file='images/pipeline-core/clean/HS-DGE_DownTid_vs_Cyte-clean.pdf', title="Human (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=False)

    #
    # MM
    #
    # Cytes vs Gonia
    #
    df = pd.read_csv('../01-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0)
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    core = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] >= 1) & (dfC['FDR_CyteGonia'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=core, file='images/pipeline-core/annotated/MM-DGE_UpCyte_vs_Gonia-ann.pdf', title="Mouse (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=True)
    plot_MA(df=df, highlight=core, file='images/pipeline-core/clean/MM-DGE_UpCyte_vs_Gonia-cln.pdf', title="Mouse (Up)Cyte vs Gonia",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=False)

    #
    # Cytes vs Tids
    #
    df = pd.read_csv('../01-diff-gene-exp/results/MM/MM-DGE_Tid_vs_Cyte.csv', index_col=0)
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    core = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] <= -1) & (dfC['FDR_TidCyte'] <= 0.05), :].index.tolist()

    plot_MA(df=df, highlight=core, file='images/pipeline-core/annotated/MM-DGE_DownTid_vs_Cyte-ann.pdf', title="Mouse (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=True)
    plot_MA(df=df, highlight=core, file='images/pipeline-core/clean/MM-DGE_DownTid_vs_Cyte-cln.pdf', title="Mouse (Down)Tid vs Cyte",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=False)

    #
    # DM
    #
    # Middle vs Apical
    #
    df = pd.read_csv('../01-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Apical.csv', index_col=0)
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    core = dfC.loc[(dfC['Middle_vs_Apical'] == True), :].index.tolist()

    plot_MA(df=df, highlight=core, file='images/pipeline-core/annotated/DM-DGE_UpMiddle_vs_Apical-ann.pdf', title="Insect (Up)Middle vs Apical",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=True)
    plot_MA(df=df, highlight=core, file='images/pipeline-core/clean/DM-DGE_UpMiddle_vs_Apical-cln.pdf', title="Insect (Up)Middle vs Apical",
            c_highlight=c_highlight, c_up='#ff9896', c_down='gray', annotate_plot=False)

    #
    # Basal vs Middle
    #
    df = pd.read_csv('../01-diff-gene-exp/results/DM/DM-DGE_Basal_vs_Middle.csv', index_col=0)
    dfC = pd.read_csv('results/pipeline-{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    core = dfC.loc[(dfC['Basal_vs_Middle'] == True), :].index.tolist()

    plot_MA(df=df, highlight=core, file='images/pipeline-core/annotated/DM-DGE_DownBasal_vs_Middle-ann.pdf', title="Insect (Down)Basal vs Middle",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=True)
    plot_MA(df=df, highlight=core, file='images/pipeline-core/clean/DM-DGE_DownBasal_vs_Middle-cln.pdf', title="Insect (Down)Basal vs Middle",
            c_highlight=c_highlight, c_up='gray', c_down='#aec7e8', annotate_plot=False)
