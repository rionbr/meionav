# coding=utf-8
# Author: Rion B Correia
# Date: Aug 06, 2019
#
# Description: Plots results of (MISS) screened DM genes
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
from matplotlib.colors import to_rgba
from pybiomart import Dataset

if __name__ == '__main__':

    df = pd.read_csv('../2-core_genes/results/DM/core_DM_meiotic_genes.csv', index_col=0)
    dfC = pd.read_csv('data/core_DM_control.csv')
    dfS = pd.read_csv('data/core_DM_screened_2019-08.csv', index_col=0, na_values='PENDING')

    # Add Gene Name
    dsDM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    dfDM = dsDM.query(attributes=['ensembl_gene_id', 'external_gene_name']).set_index('Gene stable ID')
    dfS['gene'] = dfDM['Gene name']

    # Remove NaA & PENDING genes
    dfS = dfS.loc[(~(dfS['FT1 eggs'].isnull())), :]

    # Calculations
    dfS['total-eggs'] = 0
    dfS['total-hatched'] = 0
    for ft in range(1, 5):
        col_eggs = 'FT{:d} eggs'.format(ft)
        col_hatched = 'FT{:d} hatched'.format(ft)
        col_fertate = 'FT{:d} fert-rate'.format(ft)
        dfS[col_fertate] = dfS[col_hatched] / dfS[col_eggs]
        dfS['total-eggs'] += dfS[col_eggs]
        dfS['total-hatched'] += dfS[col_hatched]

    # Mean/SD
    dfS['mean fert-rate'] = dfS[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].mean(axis=1)
    dfS['std fert-rate'] = dfS[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].std(axis=1)

    dfS = dfS.sort_values('mean fert-rate', ascending=True).reset_index()

    # Control
    dfc = dfC.groupby('control').apply(lambda x: x['hatched'] / x['eggs']).groupby('control').agg(['mean', 'std'])
    dfc = dfc.reset_index()

    miss = dfS.loc[~dfS['id_gene_DM'].isin(df.index), :]
    dfS['miss'] = ~dfS['id_gene_DM'].isin(df.index)
    dfS['color'] = dfS['miss'].map(lambda x: to_rgba('#d62728') if x==True else to_rgba('#1f77b4'))

    dfd = dfS  # All Data
    dft = dfS.iloc[:80, :]  # Data for Top axis
    dfb = dfS.iloc[80:, :]  # Data for Bottom axis

    # Plot
    fig = plt.figure(figsize=(11, 8.5))

    #
    # Three axis, two at the top and one at the bottom.
    #
    gs = gridspec.GridSpec(2, 12)
    axc = plt.subplot(gs[0, :1])  # Control
    axt = plt.subplot(gs[0, 1:])  # First part on top
    axb = plt.subplot(gs[1, :])  # Second part on bottom

    axc.errorbar(dfc.index, dfc['mean'], yerr=dfc['std'], color='#2ca02c', lw=0, elinewidth=1, capsize=2, marker='o', markersize=5)
    axt.errorbar(dft.index, dft['mean fert-rate'], yerr=dft['std fert-rate'], elinewidth=1, capsize=2, marker='o', markersize=5)
    axb.errorbar(dfb.index, dfb['mean fert-rate'], yerr=dfb['std fert-rate'], elinewidth=1, capsize=2, marker='o', markersize=5)

    dftmiss = dft.loc[ dft['miss'] == True, :]
    axt.scatter(dftmiss.index, dftmiss['mean fert-rate'], color='red', zorder=10)
    dfbmiss = dfb.loc[ dfb['miss'] == True, :]
    axb.scatter(dfbmiss.index, dfbmiss['mean fert-rate'], color='red', zorder=10)

    axc.axhline(0.75, color='#d62728', lw=1, zorder=6)
    axt.axhline(0.75, color='#d62728', lw=1, zorder=6)
    axb.axhline(0.75, color='#d62728', lw=1, zorder=6)

    axc.set_title('Control')
    axt.set_title('Screened')

    axc.set_xticks(dfc.index)
    axt.set_xticks(dft.index)
    axb.set_xticks(dfb.index)

    axc.set_xticklabels(dfc['control'], rotation=90, va='top', ha='center', fontsize='medium')
    axt.set_xticklabels(dft['gene'], rotation=90, va='top', ha='center', fontsize='medium')
    axb.set_xticklabels(dfb['gene'], rotation=90, va='top', ha='center', fontsize='x-small')

    axc.set_ylabel('Mean +/- SD Fertility Rate')
    axb.set_ylabel('Mean +/- SD Fertility Rate')
    # axt.set_xlabel('Gene')
    axb.set_xlabel('Gene')

    axc.set_xlim(dfc.index.min() - 1, dfc.index.max() + 1)
    axt.set_xlim(dft.index.min() - 1, dft.index.max() + 1)
    axb.set_xlim(dfb.index.min() - 1, dfb.index.max() + 1)
 
    axc.set_ylim(-0.02, 1.02)
    axt.set_ylim(-0.02, 1.02)
    axb.set_ylim(-0.02, 1.02)
    #
    #ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(0.1,10e5)
    axc.grid()
    axt.grid()
    axb.grid()

    # Save
    plt.tight_layout()
    fig.savefig('images/img-core_DM_miss_screened.pdf')
