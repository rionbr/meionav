# coding=utf-8
# Author: Rion B Correia
# Date: Jan 14, 2021
#
# Description: Plot the AUC curves for the ML results
#
# Instructions:
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
#
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.handlelength'] = 1.5
import matplotlib.pyplot as plt
from utils import ensurePathExists


if __name__ == '__main__':

    celltype = 'spermatocyte'
    layer = species = 'DM'

    #
    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    print('Load Embeding')
    rXYFile = 'results/unsuper-tsne/ml-unsuper-{celltype:s}-{layer:s}-tsne.csv.gz'.format(celltype=celltype, layer=layer)
    dfE = pd.read_csv(rXYFile)
    #

    print('Load X-y features')
    rMLFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-X-y.csv.gz'.format(celltype=celltype, layer=layer)
    dfX = pd.read_csv(rMLFile, index_col=0)
    # Only Core
    dfX = dfX.loc[(dfX['core'].astype(bool)), :]

    dfX['x'] = dfE['x'].tolist()
    dfX['y'] = dfE['y'].tolist()
    #print(dfX.loc[(dfX['x'] > 60) & (dfX['y'] > -255), :].sort_values('y', False).iloc[:,6:18])
    #print(dfX.loc['FBgn0010213', :])

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))

    dict_pheno_color = {
        'B': 'red',
        'A': 'blue',
        'indirect': 'lightgreen',
        'C': 'cyan',
        'D': 'magenta',
        'F': 'darkgreen',
        'G': 'violet',
        'H': 'maroon',
    }
    dfE['color'] = dfX['phenotype'].map(lambda x: dict_pheno_color[x] if x in dict_pheno_color else 'gray').to_list()
    dict_color = {
        1: 'red',
        0: 'blue'
    }
    dict_color = {
        'no-change': 'black',
        'up': 'red',
        'down': 'blue',
    }
    #dfE['color'] = dfX['mdlc-mutant-up/down'].map(lambda x: dict_color[x] if x in dict_color else 'gray').to_list()

    #
    ax.scatter(dfE['x'], dfE['y'], c=dfE['color'], s=6)

    plt.subplots_adjust(left=0.28, bottom=0.16, right=0.97, top=0.90, wspace=0.0, hspace=0.0)

    wIMGfile = 'images/unsuper-tsne/ml-{celltype:s}-{layer:s}-tsne.pdf'.format(celltype=celltype, layer=layer)
    ensurePathExists(wIMGfile)
    plt.tight_layout()
    plt.savefig(wIMGfile, dpi=300, pad_inches=0.0)

    print('Done.')
