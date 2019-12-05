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
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


if __name__ == '__main__':

    # Load genes
    dfA = pd.read_csv('../2-core_genes/results/all3-conserved/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])
    df3 = pd.read_csv('../2-core_genes/results/all3-pooling-DM/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])

    # Which genes are conserved/pooling
    df = pd.concat([dfA, df3], axis='index', join='outer').drop_duplicates()
    df['conserved'] = df.index.isin(dfA.index)
    df['pooling'] = df.index.isin(df3.index)

    # Load Screened data
    dfSc = pd.read_csv('data/conserved_DM_screened_2019-11-22.csv', index_col=0)
    dfSp = pd.read_csv('data/pooling_DM_screened_2019-11-22.csv', index_col=0)
    dfS = pd.concat([dfSc, dfSp], axis='index', join='outer').drop_duplicates()

    # Load Controls
    dfC = pd.read_csv('data/screened_DM_controls.csv')

    # Load FPKM data
    dfFPKM1 = pd.read_csv('../1-diff-gene-exp/data/DM/DM_Spermatocytes_FPKM_sample1.csv', index_col=0, usecols=['Gene ID', 'FPKM'])
    dfFPKM2 = pd.read_csv('../1-diff-gene-exp/data/DM/DM_Spermatocytes_FPKM_sample2.csv', index_col=0, usecols=['Gene ID', 'FPKM'])
    print(dfFPKM1)

    dfS_only = dfS.loc[~dfS['FT1 eggs'].isnull(), :]

    print(dfS['status'].value_counts())
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

    df['FPKM1'] = dfFPKM1['FPKM']
    df['FPKM2'] = dfFPKM2['FPKM']
    #df['FPKM1'] = df['FPKM1'].apply(np.log2)
    #df['FPKM2'] = df['FPKM2'].apply(np.log2)
    
    df = df.sort_values(['status', 'mean fert-rate'], ascending=[True, True]).reset_index()

    maxfpkm, minfpkm = df[['FPKM1','FPKM2']].max().max(), df[['FPKM1','FPKM2']].min().min()
    #print("FPKM: {:.2f}/{:.2f} (min/max)".format(minfpkm, maxfpkm))
    print(df.head())
    print(df.tail())

    # Control
    dfc = dfC.groupby('control').apply(lambda x: x['hatched'] / x['eggs']).groupby('control').agg(['mean', 'std'])
    dfc = dfc.reset_index()
    print(dfc.head())

    dfd = df  # All Data
    dft = df.iloc[:81, :]  # 81 = Data for Top axis
    df2 = df.iloc[81:361, :]  # Data for Middle axis
    df3 = df.iloc[361:641, :]  # Data for Bottom axis
    df4 = df.iloc[641:, :]  # Data for Bottom axis

    print(df.shape)
    print(dft.shape)
    print(df2.shape)
    print(df3.shape)
    print(df4.shape)
    
    # Plot
    fig = plt.figure(figsize=(11, 8.5))

    ms = 3
    """
    #
    # Two Axis
    #
    gs = gridspec.GridSpec(1, 20)
    axc = plt.subplot(gs[0, :1]) # Control
    axd = plt.subplot(gs[0, 1:]) # First part on top
    plt.setp(axd.get_yticklabels(), visible=False) # Remove yticks

    axc.errorbar(dfc.index, dfc['mean'], yerr=dfc['std'], color='#2ca02c', lw=0, elinewidth=.5, capsize=2, marker='o', ms=ms, zorder=8)
    axd.errorbar(dfd.index, dfd['mean fert-rate'], yerr=dfd['std fert-rate'], lw=0, elinewidth=.5, capsize=2, marker='o', ms=ms, zorder=8)

    axc.axhline(0.75, color='#d62728', lw=1, zorder=6)
    axd.axhline(0.75, color='#d62728', lw=1, zorder=6)

    axc.set_title('Control')
    axd.set_title('Screened')

    axc.set_xticks(dfc.index)
    #axd.set_xticks(dfd.index)

    axc.set_xticklabels(dfc['control'], rotation=90, va='top', ha='center', fontsize='medium')

    axc.set_ylabel('Mean +/- SD Fertility Rate')
    axd.set_xlabel('Gene rank')

    axc.set_xlim(dfc.index.min() - 1, dfc.index.max() + 1)
    axd.set_xlim(dfd.index.min() - 1, dfd.index.max() + 1)

    axc.set_ylim(-0.02,1.02)
    axd.set_ylim(-0.02,1.02)

    axc.grid()
    axd.grid()
    """
    #
    # Three axis, two at the top and one at the bottom.
    #
    gs = gridspec.GridSpec(14, 18, wspace=0.5, hspace=0.0)

    axc = plt.subplot(gs[0:3, :1])  # Control
    axt = plt.subplot(gs[0:3, 1:])  # First part on top
    axth = plt.subplot(gs[3, 1:])  # Heatmap below top part
    ax2 = plt.subplot(gs[5:7, :])  # Second part on bottom
    ax2h = plt.subplot(gs[7, :])
    ax3 = plt.subplot(gs[8:10, :])  # Second part on bottom
    ax3h = plt.subplot(gs[10, :])
    ax4 = plt.subplot(gs[11:13, :])  # Third part on bottom
    ax4h = plt.subplot(gs[13, :])

    for a in [axc, axt, axth, ax2, ax2h, ax3, ax3h, ax4, ax4h]:
        a.set_aspect('auto', adjustable='box', anchor='N')
        a.set_xticklabels([])
        a.set_yticklabels([])

    hnorm = mpl.colors.Normalize(vmin=minfpkm, vmax=maxfpkm)

    ebc = axc.errorbar(dfc.index, dfc['mean'], yerr=dfc['std'], color='#2ca02c', lw=0, elinewidth=1, capsize=3, marker='o', markersize=6)
    ebt = axt.errorbar(range(0, len(dft)), dft['mean fert-rate'], yerr=dft['std fert-rate'], lw=0, elinewidth=1, capsize=3, marker='o', markersize=6)
    eb2 = ax2.errorbar(range(0, len(df2)), df2['mean fert-rate'], yerr=df2['std fert-rate'], lw=0, elinewidth=1, capsize=1, marker='o', markersize=2)
    eb3 = ax3.errorbar(range(0, len(df3)), df3['mean fert-rate'], yerr=df3['std fert-rate'], lw=0, elinewidth=1, capsize=1, marker='o', markersize=2)
    eb4 = ax4.errorbar(range(0, len(df4)), df4['mean fert-rate'], yerr=df3['std fert-rate'], lw=0, elinewidth=1, capsize=1, marker='o', markersize=2)

    imth = axth.imshow(dft[['FPKM1', 'FPKM2']].T.values, cmap='Reds', norm=hnorm, aspect='equal')
    im2h = ax2h.imshow(df2[['FPKM1', 'FPKM2']].T.values, cmap='Reds', norm=hnorm, aspect='equal')
    im3h = ax3h.imshow(df3[['FPKM1', 'FPKM2']].T.values, cmap='Reds', norm=hnorm, aspect='equal')
    im4h = ax4h.imshow(df4[['FPKM1', 'FPKM2']].T.values, cmap='Reds', norm=hnorm, aspect='equal')

    axc.axhline(0.75, color='#d62728', lw=1, zorder=6)
    axt.axhline(0.75, color='#d62728', lw=1, zorder=6)
    ax2.axhline(0.75, color='#d62728', lw=1, zorder=6)
    ax3.axhline(0.75, color='#d62728', lw=1, zorder=6)

    axc.set_title('Control')
    axt.set_title('Screened')
    
    axc.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axc.set_xticks(dfc.index)
    axt.set_xticks(range(0, len(dft)))
    axth.set_yticks([0, 1])
    axth.set_xticks(range(0, len(dft)))
    ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2h.set_yticks([0, 1])
    ax2h.set_xticks(range(0, len(df2)))
    ax3.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax3h.set_yticks([0, 1])
    ax3h.set_xticks(range(0, len(df4)))
    ax4.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax4h.set_yticks([0, 1])
    ax4h.set_xticks(range(0, len(df4)))

    # ax2.set_xticks(df2.index)
    
    axc.set_xticklabels(dfc['control'], rotation=90, va='top', ha='center', fontsize='medium')
    axt.set_xticklabels([])
    axth.set_yticklabels([])
    axth.set_xticklabels(dft['gene'], rotation=90, va='top', ha='center', fontsize='medium')
    #ax2h.set_yticklabels(['1','2'])
    #ax2h.set_xticklabels([])
    #ax3h.set_xticklabels([])
    # ax2.set_xticklabels(df2['gene'], rotation=90, va='top', ha='center', fontsize='x-small')
    
    """
    plt.setp(axt.get_yticklabels(), visible=False)
    plt.setp(axt.get_xticklabels(), visible=False)
    plt.setp(ax2h.get_yticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3h.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    """
    axt.tick_params(axis='x', which='both', length=0)
    ax2.tick_params(axis='x', which='both', length=0)
    ax3.tick_params(axis='x', which='both', length=0)

    # Mark conserved in red
    for conserved, tick in zip(dft['conserved'].tolist(), axth.get_xticklabels()):
        if conserved:
            tick.set_color('red')
    axt.text(x=0.99, y=0.07, s='Genes in conserved pipeline labeled in red', transform=axt.transAxes, ha='right', va='center', zorder=20)

    axc.set_ylabel('Mean +/- SD Fertility Rate')
    ax2.set_ylabel('Cont. Mean +/- SD', fontsize='small')
    ax3.set_ylabel('Cont. Mean +/- SD', fontsize='small')
    ax4.set_ylabel('Cont. Mean +/- SD', fontsize='small')
    # axt.set_xlabel('Gene')
    ax4h.set_xlabel('Gene')

    axc.set_xlim(- 1, len(dfc))
    axt.set_xlim(- 1, len(dft))
    axth.set_xlim(-1, len(dft))
    ax2.set_xlim(-1, len(df2) + 1)
    ax2h.set_xlim(-1, len(df2) + 1)
    ax3.set_xlim(- 1, len(df3) + 1)
    ax3h.set_xlim(- 1, len(df3) + 1)
    ax4.set_xlim(- 1, len(df4) + 1)
    ax4h.set_xlim(- 1, len(df4) + 1)

    axc.set_ylim(-0.02, 1.02)
    axt.set_ylim(-0.02, 1.02)
    axth.set_ylim(-0.5, 1.5)
    ax2.set_ylim(-0.02, 1.02)
    ax2h.set_ylim(-0.5, 2)
    ax3.set_ylim(-0.02, 1.02)
    ax3h.set_ylim(-0.5, 2)
    ax4.set_ylim(-0.02, 1.02)
    ax4h.set_ylim(-0.5, 2)
    #
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(0.1,10e5)
    axc.grid()
    axt.grid()

    # Colorbar
    #cbaxes = fig.add_axes([0.01, 0.05, 0.15, 0.02])  # x, y, width, height
    #cb = plt.colorbar(imth, cax=cbaxes, orientation='horizontal')
    #cb.set_label('FPKM', rotation=0)

    # Save
    #plt.tight_layout()
    plt.subplots_adjust(left=0.06, right=0.96, bottom=0.03, top=0.94, wspace=1, hspace=0.0)
    fig.savefig('images/img-all_DM_screened.pdf')
