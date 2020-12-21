# coding=utf-8
# Author: Rion B Correia
# Date: Feb 21, 2020
#
# Description: Plots a heatmap of the pairwise jaccard comparison for conserved genes
#
# Instructions:
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
import matplotlib.gridspec as gridspec
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from itertools import combinations
from utils import ensurePathExists


if __name__ == '__main__':

    df = pd.read_csv('results/species-celltypes-threshold-stats.csv.gz', index_col=0)
    
    thresholds = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]

    # Probabilities
    df['%-genes-pc'] = df['genes-pc'] / df['genome-pc']
    df['%-genes-non-pc'] = df['genes-non-pc'] / df['genome']

    species = ['HS', 'MM', 'DM']
    #celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    celltypes = ['spermatocyte', 'enterocyte', 'neuron', 'muscle']
    species_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    
    df['species'] = pd.Categorical(df['species'], categories=species, ordered=True)
    df['celltype'] = pd.Categorical(df['celltype'], categories=celltypes, ordered=True)
    df = df.sort_values(['species', 'celltype', 'threshold'], ascending=[True, True, True])

    attr = {
        'spermatocyte': {
            'name': 'Spermatocyte',
            'color': '#d62728',
            'ls': 'solid',
        },
        'spermatogonia': {
            'name': 'Spermatogonia',
            'color': '#e377c2',
            'ls': 'solid',
        },
        'spermatid': {
            'name': 'Spermatid',
            'color': '#f7b6d2',
            'ls': 'solid',
        },
        'enterocyte': {
            'name': 'Enterocyte',
            'color': '#7f7f7f',
            'ls': 'solid',
        },
        'neuron': {
            'name': 'Neuron',
            'color': '#7f7f7f',
            'ls': 'dashed',
        },
        'muscle': {
            'name': 'Muscle',
            'color': '#7f7f7f',
            'ls': 'dotted',
        },
    }

    for pc in [True, False]:
        for specie, dfs in df.groupby('species', sort=False):
            specie_name = species_name[specie]
            print('Plot {specie:s}'.format(specie=specie_name))
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5, 4))

            # Title
            ax.set_title('{specie:s}'.format(specie=specie_name))
            handles = []
            labels = []
            for i, (celltype, dfst) in enumerate(dfs.groupby('celltype', sort=False)):

                dfst = dfs.loc[dfs['celltype'] == celltype]
                #
                color = attr[celltype]['color']
                ls = attr[celltype]['ls']
                label = attr[celltype]['name']
                #
                if specie == 'HS' and celltype == 'enterocyte':
                    label = 'Intestine'
                if specie == 'DM' and celltype == 'spermatocyte':
                    label = 'Mid testis'
                #
                if celltype in ['spermatocyte', 'spermatogonia', 'spermatid']:
                    marker = 'o'
                    ms = 6
                else:
                    marker = ''
                    ms = 0
                #
                xticks = np.arange(1, len(dfst) + 1)
                if pc:
                    handle, = ax.plot(xticks, dfst['genes-pc'].tolist(), marker=marker, ms=ms, lw=2.5, color=color, ls=ls, zorder=(12 - i))
                else:
                    handle, = ax.plot(xticks, dfst['genes-non-pc'].tolist(), marker=marker, ms=ms, lw=2.5, color=color, ls=ls, zorder=(12 - i))
                #
                ax.ticklabel_format(axis="y", style="sci", scilimits=(3, 3), useMathText=True)
                ax.set_xticks(xticks)
                ax.set_xticklabels(thresholds)
                #
                handles.append(handle)
                labels.append(label)

            if pc:
                ax.set_ylabel('Number of expressed protein coding genes')
            else:
                ax.set_ylabel('Number of expressed non-protein coding genes')
            ax.set_xlabel(r'Minimum expression threshold (TPM)')

            # Legend
            ax.legend(
                handles=handles,
                labels=labels,
                loc='upper right'
            )

            # Grid
            ax.grid()

            plt.subplots_adjust(left=0.14, right=0.97, bottom=0.12, top=0.92, wspace=0, hspace=0)
            img_path = 'images/TPM-thresholds/{specie:s}/'.format(specie=specie)
            if pc:
                file = img_path + 'img-threshold-{specie:s}-pc-genes.pdf'.format(celltype=celltype, specie=specie)
            else:
                file = img_path + 'img-threshold-{specie:s}-non-pc-genes.pdf'.format(celltype=celltype, specie=specie)
            ensurePathExists(file)
            fig.savefig(file)
            plt.close()
