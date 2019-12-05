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

    df['FPKM1'] = dfFPKM1['FPKM']
    df['FPKM2'] = dfFPKM2['FPKM']
    #df['FPKM1'] = df['FPKM1'].apply(np.log2)
    #df['FPKM2'] = df['FPKM2'].apply(np.log2)
    
    df = df.sort_values(['status', 'mean fert-rate'], ascending=[True, True]).reset_index()

    maxfpkm, minfpkm = df[['FPKM1','FPKM2']].max().max(), df[['FPKM1','FPKM2']].min().min()

    print(df.head())
    print(df.tail())
    
    # Plot
    fig = plt.figure(figsize=(11, 3.5))

    gs = gridspec.GridSpec(5, 1, wspace=0.0, hspace=0.05)
    ax = plt.subplot(gs[:3, :])
    axh = plt.subplot(gs[3:, :])

    ax.set_aspect('auto', adjustable='box', anchor='C')
    #axh.set_aspect(aspect=1, adjustable='box', anchor='N')
    axh.set_aspect(aspect=0.5, anchor='N')

    hnorm = mpl.colors.Normalize(vmin=minfpkm, vmax=maxfpkm)

    eb = ax.errorbar(range(0, len(df)), df['mean fert-rate'], yerr=df['std fert-rate'], lw=0,
        ecolor='#3182bd', elinewidth=1.0, capsize=3,
        marker='o', markersize=3.5,
        markeredgecolor='#3182bd', markeredgewidth=0.5, markerfacecolor='#6baed6', markerfacecoloralt=None)
    imh = axh.imshow(df[['FPKM1', 'FPKM2']].T.values, cmap='Reds', norm=hnorm, aspect='equal')

    ax.axhline(0.75, color='#d62728', lw=1, zorder=6)
    ax.set_title('Conserved Screened')
    
    ax.set_xticks(range(0, len(df)))
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels(np.linspace(0, 1, 5))

    
    axh.set_yticks([0, 1])
    axh.set_xticks(range(0, len(df)))
    axh.set_yticklabels(['S1','S2'], fontsize='x-small')
    axh.set_xticklabels(df['gene'], rotation=90, va='top', ha='center', fontsize='xx-small')    

    plt.setp(ax.get_xticklabels(), visible=False)

    ax.tick_params(axis='x', which='both', length=0)

    ax.set_ylabel('Mean +/- SD\nFertility Rate')
    axh.set_ylabel('FPKM\n')
    #axh.set_xlabel('Gene')

    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(-1, len(df))
    
    axh.set_ylim(-1.02, 2.02)
    axh.set_xlim(-1, len(df))    
    

    ax.grid()

    # Layout
    #plt.tight_layout()
    plt.subplots_adjust(left=0.075, right=0.99, bottom=0.07, top=0.92, wspace=0, hspace=0.0)

    # Colorbar
    cbax = fig.add_axes([0.09, 0.07, 0.12, 0.02])  # x, y, width, height
    cb = plt.colorbar(imh, cax=cbax, orientation='horizontal')
    cbax.set_title('FPKM', rotation=0, fontsize='medium')

    fig.savefig('images/img-conserved_DM_screened.pdf')
