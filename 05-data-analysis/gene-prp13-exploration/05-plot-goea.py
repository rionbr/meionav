# coding=utf-8
# Author: Rion B Correia
# Date: July 08, 2019
#
# Description: Reads GOAE results for each module and plots results
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#
from wordcloud import WordCloud
from nltk.corpus import stopwords


def plot_goea(celltype='spermatocyte', layer='DM'):

    rCSVFile = 'results/goea-{network:s}-{layer:s}.csv'.format(network=network, layer=layer)
    df = pd.read_csv(rCSVFile)

    # Trim
    df = df.loc[(df['depth'] >= 5), :]

    # All zeros are set to the smallest computable float
    df.loc[df['p_fdr_bh'] == 0.0, 'p_fdr_bh'] = np.nextafter(0, 1)
    #
    df['1-log(p)'] = 1 - (np.log(df['p_fdr_bh']))
    dft10 = df.iloc[:10, :].sort_values('1-log(p)', ascending=True)
    sl = 75 # string slice
    dft10['name'] = dft10['name'].apply(lambda x: x[0:sl] + '..' if len(x) > sl else x)

    # Plot
    fig, ax = plt.subplots(figsize=(4.7, 3.0))

    # P-values
    title = 'GOEA-{network:s}'.format(network=network)

    ind = np.arange(0, len(dft10), 1)
    bp = ax.barh(ind, 1 - np.log(dft10['p_fdr_bh']), height=0.8, facecolor=facecolor, zorder=4)
    ax.set_title(title, fontsize='large')

    minx, maxx = ax.get_xlim()
    for bar, name in zip(bp.patches, dft10['name'].tolist()):
        bx = bar.get_x()
        by = bar.get_y()
        bh = bar.get_height()
        # bw = bar.get_width()
        tx = bx + (0.01 * maxx)
        ty = (by + (bh / 2))
        ax.text(x=tx, y=ty, s=name, ha='left', va='center', fontsize='x-small', zorder=5)
    #
    #ax.axvline(x=(1 - math.log(0.01)), color='#666666', ls='dotted')
    ax.axvline(x=(1 - math.log(0.05)), color='#c7c7c7', ls='dashed')
    ax.set_yticks(ind)
    ax.set_yticklabels(dft10['GO'])
    ax.set_xlabel(r'$1 - $log($p$-value)')
    ax.set_ylim(-0.7, (10 - 0.3))
    ax.grid(axis='x', zorder=1)

    plt.subplots_adjust(left=0.21, right=0.97, bottom=0.17, top=0.89)
    #plt.tight_layout()
    #
    wIMGFile = 'images/img-goea-bars-{network:s}-{layer:s}.pdf'.format(network=network, layer=layer)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=300, bbox_inches=None, pad_inches=0.0)
    plt.close()


if __name__ == '__main__':

    network = 'rnf113'  # 'thr'
    layer = 'HS'

    dict_specie = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    facecolor = '#ffbb78'

    data_text_color = {
        1: ['ubiquitination', 'ubiquitin'],
        2: ['splicing'],
        3: ['translation', 'translational', 'cotranslational'],
        4: ['rRNA'],
        5: ['vesicle', 'transport'],
        6: ['respiration', 'respiratory', 'electron'],
        7: ['cell', 'cycle'],
        8: ['DNA', 'repair'],
        9: ['mitochondrial', 'translation', 'translational'],
        10: ['cell', 'cycle'],
        11: ['metabolic'],
        12: ['histidine', 'peptidyl', 'dephosphorylation'],
    }

    plot_goea(network, layer)
