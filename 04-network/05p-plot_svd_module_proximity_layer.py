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
mpl.rcParams['axes.labelsize'] = 'small'
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
        'Ubiquitination': 'Uniq.',
        'Translation': 'Trans.',
        'Splicing': 'Splc.',
        'rRNA processing': 'rRNA',
        'Protein catabolism': 'Prot.',
        'Mitochondrial translation': 'Mito',
        'Vesicle transport': 'Vesic.',
        'Cell division': 'Cell div.',
        'DNA repair': 'DNA rep.'
    }
    df['short-i'] = df['name-i'].replace(dict_replace)
    df['short-j'] = df['name-j'].replace(dict_replace)

    # Only plot M1 to M9
    cat_enterocyte = ['1', '2', '3', '4', '5', '6', '7', '7a', '7b', '8', '8a', '8b']
    cat_spermatocyte = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    df['id-i'] = pd.Categorical(df['id-i'].astype(str), categories=cat_enterocyte, ordered=True)
    df['id-j'] = pd.Categorical(df['id-j'].astype(str), categories=cat_spermatocyte, ordered=True)

    # Select what will be plotted
    idx_ent_hs = (df['id-i'].isin(['1', '2', '3' ,'6', '7a' ,'7b']))
    idx_ent_mm = (df['id-i'].isin(['1', '2', '3', '4', '5', '8a', '8b']))
    idx_ent_dm = (df['id-i'].isin(['2', '3', '5', '6', '7']))
    #
    idx_spe_hs = (df['id-j'].isin(['1', '2', '3', '4', '5', '6', '7', '8', '9']))
    idx_spe_mm = (df['id-j'].isin(['1', '2', '3', '4', '5', '6', '7', '8', '9']))
    idx_spe_dm = (df['id-j'].isin(['1', '2', '3', '4', '5', '6', '7', '8', '9']))

    print('Plot')
    fig = plt.figure(figsize=(3.3, 5.9))
    gs = gridspec.GridSpec(ncols=1, nrows=23, figure=fig)

    ax_hs = fig.add_subplot(gs[0:5, 0])
    ax_mm = fig.add_subplot(gs[10:16, 0])
    ax_dm = fig.add_subplot(gs[18:23, 0])

    axcb = fig.add_axes([0.04, 0.04, 0.3282, 0.009])

    title = 'Jaccard similarity on\nSVD modules across cell types'
    fig.suptitle(title)

    #
    layers = ['HS', 'MM', 'DM']
    axes = [ax_hs, ax_mm, ax_dm]
    species_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    idxslice = {
        'HS': {'enterocyte': idx_ent_hs, 'spermatocyte': idx_spe_hs},
        'MM': {'enterocyte': idx_ent_mm, 'spermatocyte': idx_spe_mm},
        'DM': {'enterocyte': idx_ent_dm, 'spermatocyte': idx_spe_dm}
    }
    index_spermatocyte = pd.MultiIndex.from_tuples([
        ('1', 'Ubiquitination'),
        ('2', 'Translation'),
        ('3', 'Splicing'),
        ('4', 'rRNA processing'),
        ('5', 'Protein catabolism'),
        ('6', 'Mitochondrial translation'),
        ('7', 'Vesicle transport'),
        ('8', 'Cell division'),
        ('9', 'DNA repair')],
        names=['id-j', 'short-j'])
    index_spermatocyte_short = pd.MultiIndex.from_tuples([
        ('1', 'Ubiq.'),
        ('2', 'Trans.'),
        ('3', 'Splc.'),
        ('4', 'rRNA'),
        ('5', 'Prot.'),
        ('6', 'Mito'),
        ('7', 'Vesic.'),
        ('8', 'Cell div.'),
        ('9', 'DNA rep.')],
        names=['id-j', 'short-j'])
    #
    cmap = cm.get_cmap('jet')
    cmap.set_bad(color='#bdbdbd')

    for layer, ax in zip(layers, axes):
        print('LAYER:', layer)

        dft = df.loc[((df['layer'] == layer) & idxslice[layer]['enterocyte'] & idxslice[layer]['spermatocyte']), :]
        print(dft)

        #
        index_enterocyte = dft[['id-i', 'name-i']].drop_duplicates().set_index(['id-i', 'name-i']).sort_index().index
        #index_spermatocyte = dft[['id-j', 'name-j']].drop_duplicates().set_index(['id-j', 'name-j']).sort_index().index

        dft = dft.pivot_table(index=['id-i', 'name-i'], columns=['id-j', 'name-j'], values='proximity', aggfunc='first')
        print(dft)
        dft = dft.reindex(index=index_enterocyte, columns=index_spermatocyte, fill_value=np.nan)
        print(dft)
        #print(dft.columns)
        #print(index_spermatocyte_short)

        im = ax.imshow(dft.values, cmap=cmap, vmin=0.0, vmax=1.0, aspect=1)

        specie = species_name[layer]
        xlabel = '{specie:s} germ cell'.format(specie=specie)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_label_position("top")

        specie = species_name[layer]
        ylabel = '{specie:s} soma'.format(specie=specie)
        ax.set_ylabel(ylabel)
        ax.yaxis.set_label_position("right")

        # Ticks
        xticks = np.arange(len(index_spermatocyte))
        yticks = np.arange(len(index_enterocyte))
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        # TickLabels
        yticklabels = ["M{cid:s}-{name:s}".format(cid=str(cid), name=name) for cid, name in index_enterocyte]
        if layer == 'HS':
            xticklabels = ["M{cid:s}-{short:s}".format(cid=str(cid), short=short) for cid, short in index_spermatocyte_short]
        else:
            xticklabels = ["M{cid:s}-{name:s}".format(cid=str(cid), name=name) for cid, name in index_spermatocyte]

        yticklabels_fake = [y.split('-')[0] for y in yticklabels]
        xticklabels_fake = [x.split('-')[0] for x in xticklabels]

        ax.set_yticklabels(yticklabels, rotation=0, fontsize='small')
        if layer == 'HS':
            ax.set_xticklabels(xticklabels, rotation=90, fontsize='small')
        else:
            ax.set_xticklabels(xticklabels_fake, rotation=90, fontsize='small')

    # Create colorbar
    cbar = plt.colorbar(im, cax=axcb, orientation='horizontal', ticks=[0, 0.5, 1], boundaries=np.linspace(0, 1, 30))
    cbar.ax.tick_params(labelsize='small')
    cbar.ax.set_title('Similarity', rotation=0, fontsize='small')

    plt.subplots_adjust(left=0.50, right=0.92, bottom=0.08, top=0.86, wspace=0.0, hspace=0.0)
    wIMGfile = 'images/proximity/img-svd-proximity-layers-{network:s}-{threshold:s}.pdf'.format(network=network, threshold=threshold_str)
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()
