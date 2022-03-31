# coding=utf-8
# Author: Rion B Correia
# Date: Jan 14, 2021
#
# Description: Plot the AUC curves for the ML results
#
# Instructions:
#
import numpy as np
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
    clftype = 'conserved-without-exp' #conserved, phenotype, conserved-without-exp
    layer = species = 'DM'

    rCSVfile = 'results/clf-{clftype:s}/{layer:s}/ml-{celltype:s}-{layer:s}-clf.csv.gz'.format(clftype=clftype, celltype=celltype, layer=layer)
    df = pd.read_csv(rCSVfile)
    #
    df['clf'] = df['clf'].replace({
        'Biased Dummy': 'Majority vote',
        'Uniform Dummy': 'Coin toss',
        'Logistic Regression': 'Logistic Reg.',
    })

    categories = ["Random Forest", "Logistic Reg.", "Linear SVM", 'Majority vote', 'Coin toss']
    df['clf'] = pd.Categorical(df['clf'], categories=categories, ordered=True)

    dfg = df.groupby('clf').agg({'mcc': ['mean', 'std']})

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4, 3))
    ind = list(range(0, len(dfg)))
    height = dfg[('mcc', 'mean')].tolist()
    bars = ax.barh(y=ind, width=height, height=0.6, zorder=4)

    facecolors = ['#d62728', '#1f77b4', '#2ca02c', '#17becf', '#bcbd22']
    edgecolors = ['#ff9896', '#aec7e8', '#98df8a', '#9edae5', '#dbdb8d']
    for bar, facecolor, edgecolor in zip(bars, facecolors, edgecolors):
        bar.set_color(facecolor)
        bar.set_edgecolor(edgecolor)
        bar.set_linewidth(1)
    ax.set_title("Matthew's corr. coef. (MCC)")
    ax.set_xlabel('score')
    ax.set_yticks(ind)
    ax.set_yticklabels(dfg.index.tolist())
    ax.invert_yaxis()
    ax.set_xlim(-0.05, 0.501)
    ax.grid()

    plt.subplots_adjust(left=0.28, bottom=0.16, right=0.97, top=0.90, wspace=0.0, hspace=0.0)

    wIMGfile = 'images/ml-auc-{clftype:s}/ml-{celltype:s}-{layer:s}-clf.pdf'.format(clftype=clftype, celltype=celltype, layer=layer)
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=300, pad_inches=0.0)

    print('Done.')
