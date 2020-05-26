# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads the similarity-layer plots the results.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from matplotlib.colors import LinearSegmentedColormap
from utils import ensurePathExists


if __name__ == '__main__':

    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    #
    print('Loading .csv')
    rCSVFile = 'results/proximity/svd-proximity-layers-{network:s}-{threshold:s}.csv.gz'.format(network=network, threshold=threshold_str)

    df = pd.read_csv(rCSVFile, index_col=0, encoding='utf-8')

    dict_replace = {
        'Ubiquitination':'Uniq.',
        'Translation': 'Tran.',
        'Splicing': 'Spli.',
        'Mitochondrial translation': 'Mito',
        'Vesicle transport': 'Vesi.',
        'rRNA processing': 'rRNA',
        'Protein catabolism': 'Prot.',
        'Cell division': 'Divi.',
        'DNA repair': 'DNA.'
    }
    df['name-i'] = df['name-i'].replace(dict_replace)
    df['name-j'] = df['name-j'].replace(dict_replace)

    # Only plot M1 to M9
    cat_enterocyte = ['1', '2', '3', '4', '5', '6', '7', '7a', '7b', '8', '8a', '8b']
    cat_spermatocyte = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    df['id-i'] = pd.Categorical(df['id-i'].astype(str), categories=cat_enterocyte, ordered=True)
    df['id-j'] = pd.Categorical(df['id-j'].astype(str), categories=cat_spermatocyte, ordered=True)

    df = df.loc[((df['id-i'] <= '8b') & (df['id-j'] <= '9')), :]

    print('Plot')
    fig = plt.figure(figsize=(8, 3.0))
    gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)
    ax_hs = fig.add_subplot(gs[0, 0])
    ax_mm = fig.add_subplot(gs[0, 1])
    ax_dm = fig.add_subplot(gs[0, 2])

    axcb = fig.add_axes([0.06, 0.12, 0.14, 0.021])

    title = 'Jaccard similarity on SVD modules across cell types'
    fig.suptitle(title)

    #
    layers = ['HS', 'MM', 'DM']
    axes = [ax_hs, ax_mm, ax_dm]
    species_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    #
    cmap = cm.get_cmap('jet')
    cmap.set_bad(color='#bdbdbd')

    for layer, ax in zip(layers, axes):
        print('LAYER:', layer)
        #
        dft = df.loc[(df['layer'] == layer), :]
        #
        index_enterocyte = dft[['id-i', 'name-i']].drop_duplicates().set_index(['id-i', 'name-i']).sort_index().index
        index_spermatocyte = dft[['id-j', 'name-j']].drop_duplicates().set_index(['id-j', 'name-j']).sort_index().index
        #
        dft = dft.pivot_table(index=['id-i', 'name-i'], columns=['id-j', 'name-j'], values='proximity', aggfunc='first')
        dft = dft.reindex(index=index_enterocyte, columns=index_spermatocyte, fill_value=np.nan)

        im = ax.imshow(dft.values, cmap=cmap, vmin=0.0, vmax=1.0)

        specie = species_name[layer]
        xlabel = '{specie:s} spermatocyte'.format(specie=specie)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_label_position("top")

        specie = species_name[layer]
        ylabel = '{specie:s} enterocyte'.format(specie=specie)
        ax.set_ylabel(ylabel)
        ax.yaxis.set_label_position("right")

        # Ticks
        xticks = np.arange(len(index_spermatocyte))
        yticks = np.arange(len(index_enterocyte))
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        # TickLabels
        yticklabels = ["M{cid:s}-{name:s}".format(cid=str(cid), name=name) for cid, name in index_enterocyte]
        xticklabels = ["M{cid:s}-{name:s}".format(cid=str(cid), name=name) for cid, name in index_spermatocyte]
        yticklabels_fake = [y.split('-')[0] for y in yticklabels]
        xticklabels_fake = [x.split('-')[0] for x in xticklabels]

        ax.set_yticklabels(yticklabels, rotation=0, fontsize='small')
        ax.set_xticklabels(xticklabels, rotation=90, fontsize='small')

    # Create colorbar
    cbar = plt.colorbar(im, cax=axcb, orientation='horizontal', ticks=[0, 0.5, 1], boundaries=np.linspace(0, 1, 30))
    cbar.ax.tick_params(labelsize='small')
    cbar.ax.set_title('Similarity', rotation=0, fontsize='small')

    plt.subplots_adjust(left=0.1, right=0.96, bottom=0.25, top=0.98, wspace=0.55, hspace=0.0)
    wIMGfile = 'images/proximity/img-svd-proximity-layers-{network:s}-{threshold:s}.pdf'.format(network=network, threshold=threshold_str)
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()
