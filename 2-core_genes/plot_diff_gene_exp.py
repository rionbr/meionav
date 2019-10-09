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
            c_up='black', c_up_core='black', c_up_pool='black',
            c_not='black',
            c_down='black', c_down_core='black', c_down_pool='black',
            m_up='o', m_up_core='P', m_up_pool='X',
            m_not='o',
            m_down='o', m_down_core='P', m_down_pool='X',
            ):
    s = 5
    lw = 0
    alpha = 0.8
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

    minLogFC = math.log2(2)

    # Mask Data
    mask_core = (df.index.isin(core))
    mask_pool = (df.index.isin(pool))
    mask_logfc = (df['logFC'].abs() >= minLogFC)
    mask_fdr = (df['FDR'] <= 0.05)
    #
    mask_up_core = mask_core & mask_logfc & mask_fdr & (df['logFC'] >= 0)
    mask_up_pool = mask_pool & ~mask_core & mask_logfc & mask_fdr & (df['logFC'] >= 0)
    mask_up = ~mask_core & ~mask_pool & mask_logfc & mask_fdr & (df['logFC'] >= 0)
    mask_down_core = mask_core & mask_logfc & mask_fdr & (df['logFC'] <= 0)
    mask_down_pool = mask_pool & ~mask_core & mask_logfc & mask_fdr & (df['logFC'] <= 0)
    mask_down = ~mask_core & ~mask_pool & mask_logfc & mask_fdr & (df['logFC'] <= 0)
    mask_not = ~(mask_up | mask_up_core | mask_up_pool | mask_down | mask_down_core | mask_down_pool)
    # Filter Data
    dfUp = df.loc[mask_up, :]
    dfUpC = df.loc[mask_up_core, :]
    dfUpP = df.loc[mask_up_pool, :]
    dfNot = df.loc[mask_not, :]
    dfDown = df.loc[mask_down, :]
    dfDownC = df.loc[mask_down_core, :]
    dfDownP = df.loc[mask_down_pool, :]

    # Sanity Check
    if not set(df.index.tolist()) == set(dfUp.index.tolist() + dfUpC.index.tolist() + dfUpP.index.tolist() + dfNot.index.tolist() + dfDown.index.tolist() + dfDownC.index.tolist() + dfDownP.index.tolist()):
        raise ValueError('Sanity check has failed!')
    # Counts
    n_up, n_up_core, n_up_pool, n_not, n_down, n_down_core, n_down_pool = len(dfUp), len(dfUpC), len(dfUpP), len(dfNot), len(dfDown), len(dfDownC), len(dfDownP)
    # Plot
    ax.scatter(dfUp['logCPM'], dfUp['logFC'], c=c_up, s=s, lw=lw, alpha=alpha, marker=m_up, zorder=2)
    ax.scatter(dfUpC['logCPM'], dfUpC['logFC'], c=c_up_core, s=s * 2, lw=lw, alpha=alpha, marker=m_up_core, zorder=4)
    ax.scatter(dfUpP['logCPM'], dfUpP['logFC'], c=c_up_pool, s=s * 2, lw=lw, alpha=alpha, marker=m_up_pool, zorder=3)
    ax.scatter(dfNot['logCPM'], dfNot['logFC'], c=c_not, s=s / 3, lw=lw, alpha=alpha, marker=m_not, zorder=1)
    ax.scatter(dfDown['logCPM'], dfDown['logFC'], c=c_down, s=s, lw=lw, alpha=alpha, marker=m_down, zorder=2)
    ax.scatter(dfDownC['logCPM'], dfDownC['logFC'], c=c_down_core, s=s * 2, lw=lw, alpha=alpha, marker=m_down_core, zorder=4)
    ax.scatter(dfDownP['logCPM'], dfDownP['logFC'], c=c_down_pool, s=s * 2, lw=lw, alpha=alpha, marker=m_down_pool, zorder=3)
    # Draw a line at y=(-1,0,1)
    ax.axhline(y=-1, color='b', lw=1, linestyle='--')
    ax.axhline(y=0, color='gray', lw=1, linestyle='--')
    ax.axhline(y=+1, color='b', lw=1, linestyle='--')

    ax.set_xlim(-1,18)
    ax.set_ylim(-15,15)

    # Number of Selected Genes
    # Up
    if n_up_core > 0:
        strings = [
            '{:,d}'.format(n_up_core),
            '+',
            '{:,d}'.format(n_up_pool),
            '+',
            '{:,d}'.format(n_up),
            '=',
            '{:,d}'.format(n_up_core + n_up_pool + n_up)]
        colors = [c_up_core, 'black', c_up_pool, 'black', c_up, 'black', 'black']
    else:
        strings = ['{:,d}'.format(n_up)]
        colors = ['black']
    rainbow_text(ax=ax, x=8, y=13, strings=strings, colors=colors, ha='left', va='center', fontsize='large')
    # Not
    ax.text(x=14, y=-0.1, s='{:,d}'.format(n_not), color=c_not, ha='left', va='center', fontsize='large')
    # Down
    if n_down_core > 0:
        strings = [
            '{:,d}'.format(n_down_core),
            '+',
            '{:,d}'.format(n_down_pool),
            '+',
            '{:,d}'.format(n_down),
            '=',
            '{:,d}'.format(n_down_core + n_down_pool + n_down)]
        colors = [c_down_core, 'black', c_down_pool, 'black', c_down, 'black', 'black']
    else:
        strings = ['{:,d}'.format(n_down)]
        colors = ['black']
    rainbow_text(ax=ax, x=8, y=-13, strings=strings, colors=colors, ha='left', va='center', fontsize='large')
    
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

    core_pipeline = 'all3-conserved-FDRp05'
    pool_pipeline = 'all3-pooling-DM-FDRp01'
    #
    # [H]omo [S]apiens
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0, nrows=None)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] > 0), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0, nrows=None)
    pool = dfP.loc[(dfP['Cyte_vs_Gonia']== True) & (dfP['logFC_CyteGonia'] > 0), :].index.tolist()
    print(len(core))
    print(len(pool))
    plot_MA(df=df, core=core, pool=pool, file='images/HS-DGE_UpCyte_vs_Gonia.pdf', title="HS (Up)Cyte vs Gonia",
        x_up=10, y_up=5, x_not=15, y_not=0, x_down=15, y_down=-7,
        c_up_core='#d62728', c_up_pool='#9467bd', c_up='#ff9896', c_down='gray'
    )
    #
    # Cytes vs Tids (interested in genes downregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Tid_vs_Cyte.csv', index_col=0, nrows=None)
    df.index = df.index.map(lambda x: x.split('.')[0])
    dfC = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] < 0), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0, nrows=None)
    pool = dfP.loc[(dfP['Tid_vs_Cyte'] == True) & (dfP['logFC_TidCyte'] < 0), :].index.tolist()
    
    plot_MA(df=df, core=core, pool=pool, file='images/HS-DGE_DownTid_vs_Cyte.pdf', title="HS (Down)Tid vs Cyte",
        x_up=14, y_up=6, x_not=14, y_not=0, x_down=9, y_down=-13,
        c_down_core='#1f77b4', c_down_pool='#9467bd', c_up='gray', c_down='#aec7e8'
    )
    #
    # MM
    #
    # Cytes vs Gonia (interested in genes upregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['Cyte_vs_Gonia'] == True) & (dfC['logFC_CyteGonia'] > 0), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0, nrows=None)
    pool = dfP.loc[(dfP['Cyte_vs_Gonia']== True) & (dfP['logFC_CyteGonia'] > 0), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/MM-DGE_UpCyte_vs_Gonia.pdf', title="MM (Up)Cyte vs Gonia",
        x_up=1, y_up=9, x_not=11, y_not=0, x_down=1, y_down=-9,
        c_up_core='#d62728', c_up_pool='#9467bd', c_up='#ff9896', c_down='gray'
    )
    #
    # Cytes vs Tids (interested in genes downregulated in Cytes)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Tid_vs_Cyte.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['Tid_vs_Cyte'] == True) & (dfC['logFC_TidCyte'] < 0), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0, nrows=None)
    pool = dfP.loc[(dfP['Tid_vs_Cyte'] == True) & (dfP['logFC_TidCyte'] < 0), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/MM-DGE_DownTid_vs_Cyte.pdf', title="MM (Down)Tid vs Cyte",
        x_up=11, y_up=4, x_not=11, y_not=0, x_down=7, y_down=-6,
        c_down_core='#1f77b4', c_down_pool='#9467bd', c_down='#aec7e8', c_up='gray'
    )
    #
    # DM
    #
    # Middle vs Apical (interested in genes upregulated in Middle)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Apical.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['Middle_vs_Apical'] == True) & (dfC['logFC_MiddleApical'] > 0), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0, nrows=None)
    pool = dfP.loc[(dfP['Middle_vs_Apical']== True) & (dfP['logFC_MiddleApical'] > 0), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/DM-DGE_UpMiddle_vs_Apical.pdf', title="DM (Up)Middle vs Apical",
        x_up=8, y_up=4, x_not=12, y_not=-0.40, x_down=12, y_down=-5,
        c_up_core='#d62728', c_up_pool='#9467bd', c_up='#ff9896', c_down='gray'
    )
    #
    # Middle vs Basal (interested in genes downregulated in Middle)
    #
    df = pd.read_csv('../1-diff-gene-exp/results/DM/DM-DGE_Basal_vs_Middle.csv', index_col=0, nrows=None)
    dfC = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=core_pipeline), index_col=0, nrows=None)
    core = dfC.loc[(dfC['Basal_vs_Middle'] == True) & (dfC['logFC_BasalMiddle'] < 0), :].index.tolist()
    dfP = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pool_pipeline), index_col=0, nrows=None)
    pool = dfP.loc[(dfP['Basal_vs_Middle'] == True) & (dfP['logFC_BasalMiddle'] < 0), :].index.tolist()

    plot_MA(df=df, core=core, pool=pool, file='images/DM-DGE_DownBasal_vs_Middle.pdf', title="DM (Down)Basal vs Middle",
        x_up=12.5, y_up=4, x_not=12, y_not=0.375, x_down=9, y_down=-6,
        c_down_core='#1f77b4', c_down='#aec7e8', c_up='gray'
    )
