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
    dfA = pd.read_csv('../2-core_genes/results/all3-conserved-FDRp05/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])
    dfB = pd.read_csv('../2-core_genes/results/all3-pooling-DM-FDRp01/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])
    
    # Which genes are conserved/pooling
    df = pd.concat([dfA, dfB], axis='index', join='outer').drop_duplicates()
    df['conserved'] = df.index.isin(dfA.index)
    df['pooling'] = df.index.isin(dfB.index)

    # Load Screened data
    dfS = pd.read_csv('data/core_DM_screened_2019-09.csv', index_col=0, na_values='PENDING')
    dfC = pd.read_csv('data/core_DM_control.csv')

    print(dfS.head())
    list_pending = dfS.loc[dfS['FlyBase Genotype'] == 'PENDING', []].index
    dfSonly = dfS.loc[~dfS['FT1 eggs'].isnull(), :]
    list_screened = dfSonly.index

    def map_status(x, list_pending, list_screened):
        if x.name in list_pending:
            return "SCHEDULED"
        elif x.name in list_screened:
            return "SCREENED"
        else:
            return "PENDING"

    df['status'] = df.apply(map_status, args=(list_pending, list_screened), axis='columns')
    df['status'] = pd.Categorical(df['status'], categories=['SCREENED', 'SCHEDULED', 'PENDING'], ordered=True)

    cols = ['FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched']
    df[cols] = dfSonly[cols]
    
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

    df = df.sort_values(['status', 'mean fert-rate'], ascending=[True,True]).reset_index()

    # Control
    dfc = dfC.groupby('control').apply(lambda x: x['hatched'] / x['eggs']).groupby('control').agg(['mean', 'std'])
    dfc = dfc.reset_index()
    print(dfc.head())

    dfd = df # All Data
    dft = df.iloc[:70, :] # Data for Top axis
    dfb = df.iloc[70:, :] # Data for Bottom axis

    # Plot
    fig = plt.figure(figsize=(11, 7))

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
    gs = gridspec.GridSpec(2, 12)
    axc = plt.subplot(gs[0, :1]) # Control
    axt = plt.subplot(gs[0, 1:]) # First part on top
    axb = plt.subplot(gs[1, :]) # Second part on bottom

    axc.errorbar(dfc.index, dfc['mean'], yerr=dfc['std'], color='#2ca02c', lw=0, elinewidth=1, capsize=3, marker='o', markersize=6)
    axt.errorbar(dft.index, dft['mean fert-rate'], yerr=dft['std fert-rate'], lw=0, elinewidth=1, capsize=3, marker='o', markersize=6)
    axb.errorbar(dfb.index, dfb['mean fert-rate'], yerr=dfb['std fert-rate'], lw=0, elinewidth=1, capsize=1, marker='o', markersize=2)

    axc.axhline(0.75, color='#d62728', lw=1, zorder=6)
    axt.axhline(0.75, color='#d62728', lw=1, zorder=6)
    axb.axhline(0.75, color='#d62728', lw=1, zorder=6)

    axc.set_title('Control')
    axt.set_title('Screened')

    axc.set_xticks(dfc.index)
    axt.set_xticks(dft.index)
    #axb.set_xticks(dfb.index)

    axc.set_xticklabels(dfc['control'], rotation=90, va='top', ha='center', fontsize='medium')
    axt.set_xticklabels(dft['gene'], rotation=90, va='top', ha='center', fontsize='medium')    
    #axb.set_xticklabels(dfb['gene'], rotation=90, va='top', ha='center', fontsize='x-small')
    
    # Mark conserved in red
    for conserved, tick in zip(dft['conserved'].tolist(), axt.get_xticklabels()):
        if conserved:
            tick.set_color('red')
    axt.text(x=0.99, y=0.07, s='Genes in conserved pipeline labeled in red', transform=axt.transAxes, ha='right', va='center', zorder=20)

    axc.set_ylabel('Mean +/- SD Fertility Rate')
    axb.set_ylabel('Mean +/- SD Fertility Rate')
    #axt.set_xlabel('Gene')
    axb.set_xlabel('Gene')

    axc.set_xlim(dfc.index.min() - 1, dfc.index.max() + 1)
    axt.set_xlim(dft.index.min() - 1, dft.index.max() + 1)
    axb.set_xlim(dfb.index.min() - 1, dfb.index.max() + 1)
 
    axc.set_ylim(-0.02,1.02)
    axt.set_ylim(-0.02,1.02)
    axb.set_ylim(-0.02,1.02)
    #
    #ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(0.1,10e5)
    axc.grid()
    axt.grid()
    axb.grid()

    # Save
    plt.tight_layout()
    fig.savefig('images/img-all_DM_screened.pdf')
