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

    # Keep only genes in core list
    dfS = dfS.loc[dfS.index.isin(df.index), :]
    print(dfS.iloc[:10,:10])

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

    # Plot
    fig = plt.figure(figsize=(11, 3.5))
    ms = 6
    elinewidth = 1
    capsize = 3
    #
    # Two Axis
    #
    gs = gridspec.GridSpec(1, 10)
    axc = plt.subplot(gs[0, :1])  # Control
    axd = plt.subplot(gs[0, 1:])  # First part on top
    plt.setp(axd.get_yticklabels(), visible=False)  # Remove yticks

    axc.errorbar(dfc.index, dfc['mean'], yerr=dfc['std'], color='#2ca02c', lw=0, elinewidth=elinewidth, capsize=capsize, marker='o', ms=ms, zorder=8)
    axd.errorbar(dfd.index, dfd['mean fert-rate'], yerr=dfd['std fert-rate'], lw=0, elinewidth=elinewidth, capsize=capsize, marker='o', ms=ms, zorder=8)

    axc.axhline(0.75, color='#d62728', lw=1, zorder=6)
    axd.axhline(0.75, color='#d62728', lw=1, zorder=6)

    axc.set_title('Control')
    axd.set_title('Screened')

    axc.set_xticks(dfc.index)
    axd.set_xticks(dfd.index)

    axc.set_xticklabels(dfc['control'], rotation=90, va='top', ha='center', fontsize='medium')
    axd.set_xticklabels(dfd['gene'], rotation=90, va='top', ha='center', fontsize='medium')

    axc.set_ylabel('Mean +/- SD Fertility Rate')

    axc.set_xlim(dfc.index.min() - 1, dfc.index.max() + 1)
    axd.set_xlim(dfd.index.min() - 1, dfd.index.max() + 1)
 
    axc.set_ylim(-0.02,1.02)
    axd.set_ylim(-0.02,1.02)

    axc.grid()
    axd.grid()

    # Save
    plt.tight_layout()
    fig.savefig('images/img-core_DM_screened.pdf')
