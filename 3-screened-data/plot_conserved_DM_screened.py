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
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


if __name__ == '__main__':

    # Load genes
    df = pd.read_csv('../2-core_genes/results/all3-conserved/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])

    # Load Screened data
    dfSc = pd.read_csv('data/conserved_DM_screened_2019-11-22.csv', index_col=0)
    dfSp = pd.read_csv('data/pooling_DM_screened_2019-11-22.csv', index_col=0)
    dfS = pd.concat([dfSc, dfSp], axis='index', join='outer').drop_duplicates()

    # Load FPKM data
    dfFPKM1 = pd.read_csv('../1-diff-gene-exp/data/DM/DM_Spermatocytes_FPKM_sample1.csv', index_col=0, usecols=['Gene ID', 'FPKM'])
    dfFPKM2 = pd.read_csv('../1-diff-gene-exp/data/DM/DM_Spermatocytes_FPKM_sample2.csv', index_col=0, usecols=['Gene ID', 'FPKM'])

    dfS_only = dfS.loc[~dfS['FT1 eggs'].isnull(), :]

    status_cats = ['Screened', 'Crossed 21/11', 'Crossed 11/11', 'Order v27870', 'Ordered', 'Pending']
    dfS['status'] = pd.Categorical(dfS['status'], categories=status_cats, ordered=True)
    df['status'] = dfS['status']

    cols = ['FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched']
    df[cols] = dfS_only[cols]

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

    # Previously known to affect infertility in DS?
    df['prev'] = dfS['Previously known to affect male fertility/sperm cells'].apply(lambda x: 1 if 'Yes' in x else 0)
    df['RNAi'] = dfS['Previous ref to RNAi working'].apply(lambda x: 1 if 'Yes' in x else 0)

    # FPKM
    df['FPKM1'] = dfFPKM1['FPKM']
    df['FPKM2'] = dfFPKM2['FPKM']
    df['FPKM'] = df[['FPKM1', 'FPKM2']].mean(axis='columns')
    df['logFPKM'] = df['FPKM'].apply(lambda x: np.log2(x + 1))

    df = df.sort_values(['status', 'mean fert-rate'], ascending=[False, False]).reset_index()

    maxfpkm, minfpkm = df['logFPKM'].max(), df['logFPKM'].min()

    print(df.head())
    print(df.tail())
    print("logFPKM: {:.2f}/{:.2f}".format(minfpkm, maxfpkm))
    # Plot
    fig = plt.figure(figsize=(3.5, 11))

    gs = gridspec.GridSpec(nrows=1, ncols=15, wspace=0.5, hspace=0.0)
    
    axr = plt.subplot(gs[:, 0:5])
    axp = plt.subplot(gs[:, 5])
    axh = plt.subplot(gs[:, 6])
    ax = plt.subplot(gs[:, 7:14])

    axr.set_aspect(aspect='auto', anchor='E')
    axp.set_aspect(aspect='auto', anchor='E')
    axh.set_aspect(aspect='auto', anchor='E')
    ax.set_aspect(aspect='auto', adjustable='box', anchor='W')
    
    hnorm = mpl.colors.Normalize(vmin=minfpkm, vmax=maxfpkm)
    pcmap = colors.ListedColormap(['white', '#9467bd'])
    pbounds = [0, 0.5, 1]
    pnorm = colors.BoundaryNorm(pbounds, pcmap.N)
    rcmap = colors.ListedColormap(['white', '#17becf'])
    rbounds = [0, 0.5, 1]
    rnorm = colors.BoundaryNorm(rbounds, rcmap.N)

    # Plot
    eb = ax.errorbar(df['mean fert-rate'], range(0, len(df)), xerr=df['std fert-rate'], lw=0,
                     ecolor='#3182bd', elinewidth=1.0, capsize=3,
                     marker='o', markersize=3.5,
                     markeredgecolor='#3182bd', markeredgewidth=0.5,
                     markerfacecolor='#6baed6', markerfacecoloralt=None
                     )

    imr = axr.imshow(df[['RNAi']].values, cmap=rcmap, norm=rnorm, aspect='equal')
    imp = axp.imshow(df[['prev']].values, cmap=pcmap, norm=pnorm)
    imh = axh.imshow(df[['logFPKM']].values, cmap='Reds', norm=hnorm, aspect='equal')

    # Horizontal Line
    ax.axvline(0.75, color='#d62728', lw=1, zorder=6)

    # axr
    axr.set_xticks([0])
    axr.set_yticks(range(0, len(df)))
    axr.set_xticklabels(['RNAi efficiency'], fontsize='small', rotation=35, ha='right')
    axr.set_yticklabels(df['gene'], rotation=0, va='center', ha='right', fontsize='xx-small')
    axr.set_xlim(-1.02, 1.02)
    axr.set_ylim(-1, len(df))
    # axp
    axp.set_xticks([0])
    axp.set_yticks([])
    axp.set_xticklabels(['Prev. reported'], fontsize='small', rotation=35, ha='right')
    axp.set_yticklabels([])
    plt.setp(axp.get_yticklabels(), visible=False)
    axp.set_xlim(-1.02, 1.02)
    axp.set_ylim(-1, len(df))

    # axh
    axh.set_xticks([0])
    axh.set_yticks([])
    axh.set_xticklabels(['FPKM'], fontsize='small', rotation=35, ha='right')
    axp.set_yticklabels([])
    plt.setp(axh.get_yticklabels(), visible=False)
    axh.set_xlim(-1.02, 1.02)
    axh.set_ylim(-1, len(df))

    # ax
    fig.suptitle('Core metazoan meiotic genes')
    ax.set_xlabel('Fertility Rate (Mean +/- SD)', fontsize='small')
    ax.set_yticks(range(0, len(df)))
    ax.set_xticks(np.linspace(0, 1, 5))
    ax.set_xticklabels(np.linspace(0, 1, 5), fontsize='small', rotation=0)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='y', which='both', length=0)
    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-1, len(df))
    ax.grid(linewidth=0.5)

    # Layout
    plt.subplots_adjust(left=0.09, right=0.98, bottom=0.07, top=0.96, wspace=0, hspace=0.0)

    # RNAi Colorbar
    rnai_cbax = fig.add_axes([0.02, 0.07, 0.02, 0.06])  # x, y, width, height
    rnai_cb = plt.colorbar(imr, cax=rnai_cbax, ticks=[-1, -0.5, 0, 0.5, 1], orientation='vertical')
    rnai_cbax.set_title(r'RNAi efficiency', rotation=90, fontsize='small')
    rnai_cbax.yaxis.set_label_position('left')
    rnai_cb.set_ticklabels(['No', '', 'Yes'])
    rnai_cb.ax.tick_params(labelsize='small')
    cb_ticks = rnai_cb.ax.get_yticklabels()
    print(cb_ticks)
    cb_ticks[0].set_verticalalignment('bottom')
    cb_ticks[-1].set_verticalalignment('top')

    # Prev Colorbar
    prev_cbax = fig.add_axes([0.02, 0.24, 0.02, 0.06])  # x, y, width, height
    prev_cb = plt.colorbar(imp, cax=prev_cbax, ticks=[-1, -0.5, 0, 0.5, 1], orientation='vertical')
    prev_cbax.set_title(r'Previously reported', rotation=90, fontsize='small')
    prev_cbax.yaxis.set_label_position('left')
    prev_cb.set_ticklabels(['No', '', 'Yes'])
    prev_cb.ax.tick_params(labelsize='small')
    cb_ticks = prev_cb.ax.get_yticklabels()
    cb_ticks[0].set_verticalalignment('bottom')
    cb_ticks[-1].set_verticalalignment('top')

    # FPKM Colorbar
    fpkm_cbax = fig.add_axes([0.02, 0.44, 0.02, 0.12])  # x, y, width, height
    fpkm_cb = plt.colorbar(imh, cax=fpkm_cbax, orientation='vertical')
    fpkm_cbax.set_title(r'log$_2$(FPKM+1)', rotation=90, fontsize='small')
    fpkm_cbax.yaxis.set_label_position('left')
    fpkm_cb.ax.tick_params(labelsize='small')

    fig.savefig('images/img-conserved_DM_screened.pdf')
