# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads the similarity-celltype and plots the results.
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
from utils import ensurePathExists


if __name__ == '__main__':

    celltype = 'spermatocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    dict_replace = {
        'Ubiquitination': 'Ubiq.',
        'Splicing': 'Splc.',
        'Translation': 'Trsl.',
        'rRNA regulation': 'rRNA',
        'Vesicle transport': 'Ves. trsp.',
        'Respiration': 'Resp.',
        'Cell cycle': 'Cell cyc.',
        'DNA repair': 'DNA rep.',
        'Mitochondrial translation': 'M. trsl.',
        'Cell cycle (II)': 'Cell cyc. II',
        'Metabolism': 'Metb.',
        'Peptidyl-histidine dephosphorylation': 'Pep. dephospho.'
    }
    #
    layers = ['HS', 'MM', 'DM']
    layer_pairs = [('HS', 'MM'), ('HS', 'DM'), ('MM', 'DM')]
    species_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    #
    fig, axes = plt.subplots(figsize=(8.5, 7.2), nrows=3, ncols=3, sharex=True, sharey=True)
    #
    cmap = cm.get_cmap('jet')
    cmap.set_bad(color='#bdbdbd')

    #
    # First Column
    #
    print('Loading spermatocyte csv')
    celltype = 'spermatocyte'
    celltype_str = 'meiotic cell' if celltype == 'spermatocyte' else 'soma'
    rCSVFile = 'results/module-proximity/module-proximity-{celltype:s}-{network:s}-{threshold:s}.csv.gz'.format(celltype='spermatocyte', network=network, threshold=threshold_str)
    df = pd.read_csv(rCSVFile, index_col=0)
    # Only plot M1 - M6
    df = df.loc[((df['id-i'] <= 6) & (df['id-j'] <= 6)), :]
    #
    df['short-i'] = df['name-i'].replace(dict_replace)
    df['short-j'] = df['name-j'].replace(dict_replace)
    #
    for (layer_i, layer_j), ax in zip(layer_pairs, [axes[0][0], axes[1][0], axes[2][0]]):
        #
        dft = df.loc[(df['layer-i'] == layer_i) & (df['layer-j'] == layer_j), :]

        dft = dft.pivot_table(index=['id-i', 'name-i', 'short-i'], columns=['id-j', 'name-j', 'short-j'], values='proximity', aggfunc='first')

        index = dft.index if len(dft.index) > len(dft.columns) else dft.columns
        dft = dft.reindex(index=index, columns=index, fill_value=np.nan).T

        xticks = yticks = list(range(max(dft.shape)))
        yticklabels = ['M{id:d}-{name:s}'.format(id=id, name=name) for id, name, short in dft.index.values]
        xticklabels = ['M{id:d}'.format(id=id) for id, name, short in dft.index.values]

        im = ax.imshow(dft.values, cmap=cmap, vmin=0.0, vmax=1.0)

        species_i = species_name[layer_i]
        xlabel = '{species:s} {celltype:s}'.format(species=species_i, celltype=celltype_str)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_label_position("top")

        species_j = species_name[layer_j]
        ylabel = '{species:s} {celltype:s}'.format(species=species_j, celltype=celltype_str)
        ax.set_ylabel(ylabel)
        ax.yaxis.set_label_position("right")

        # Ticks
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        ax.set_yticklabels(yticklabels, rotation=0, ha='left')
        ax.get_yaxis().set_tick_params(pad=100)
    axes[2][0].set_xticklabels(xticklabels, rotation=0)

    #
    # Second Column
    #
    print('Loading enterocyte csv')
    celltype = 'enterocyte'
    celltype_str = 'meiotic cell' if celltype == 'spermatocyte' else 'soma'
    rCSVFile = 'results/module-proximity/module-proximity-{celltype:s}-{network:s}-{threshold:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str)
    df = pd.read_csv(rCSVFile, index_col=0)
    # Only plot M1 - M6
    df = df.loc[((df['id-i'] <= 6) & (df['id-j'] <= 6)), :]
    #
    df['short-i'] = df['name-i'].replace(dict_replace)
    df['short-j'] = df['name-j'].replace(dict_replace)
    #
    for (layer_i, layer_j), ax in zip(layer_pairs, [axes[0][1], axes[1][1], axes[2][1]]):
        #
        dft = df.loc[(df['layer-i'] == layer_i) & (df['layer-j'] == layer_j), :]

        dft = dft.pivot_table(index=['id-i', 'name-i', 'short-i'], columns=['id-j', 'name-j', 'short-j'], values='proximity', aggfunc='first')

        index = dft.index if len(dft.index) > len(dft.columns) else dft.columns
        dft = dft.reindex(index=index, columns=index, fill_value=np.nan).T

        xticks = yticks = list(range(max(dft.shape)))
        yticklabels = ['M{id:d}-{name:s}'.format(id=id, name=name) for id, name, short in dft.index.values]
        xticklabels = ['M{id:d}'.format(id=id) for id, name, short in dft.index.values]

        im = ax.imshow(dft.values, cmap=cmap, vmin=0.0, vmax=1.0)

        species_i = species_name[layer_i]
        xlabel = '{species:s} {celltype:s}'.format(species=species_i, celltype=celltype_str)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_label_position("top")

        species_j = species_name[layer_j]
        ylabel = '{species:s} {celltype:s}'.format(species=species_j, celltype=celltype_str)
        ax.set_ylabel(ylabel)
        ax.yaxis.set_label_position("right")

        # Ticks
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        ax.set_yticklabels(yticklabels, rotation=0)
        ax.set_xticklabels(xticklabels, rotation=0)


    #
    # Third Column
    #
    rCSVFile = 'results/module-proximity/module-proximity-layers-{network:s}-{threshold:s}.csv.gz'.format(network=network, threshold=threshold_str)
    df = pd.read_csv(rCSVFile, index_col=0, encoding='utf-8')
    # Only plot M1 - M6
    df = df.loc[((df['id-i'] <= 6) & (df['id-j'] <= 6)), :]
    #
    df['short-i'] = df['name-i'].replace(dict_replace)
    df['short-j'] = df['name-j'].replace(dict_replace)
    #
    for layer, ax in zip(layers, [axes[0][2], axes[1][2], axes[2][2]]):
        print('Layer:', layer)

        dft = df.loc[(
            (df['layer'] == layer)
        ), :]

        # Pivot
        dft = dft.pivot_table(index=['id-i', 'name-i', 'short-i'], columns=['id-j', 'name-j', 'short-j'], values='proximity', aggfunc='first')

        index = dft.index if len(dft.index) > len(dft.columns) else dft.columns
        dft = dft.reindex(index=index, columns=index, fill_value=np.nan)

        xticks = yticks = list(range(max(dft.shape)))
        yticklabels = ['M{id:d}-{name:s}'.format(id=id, name=name) for id, name, short in dft.index.values]
        xticklabels = ['M{id:d}'.format(id=id) for id, name, short in dft.index.values]

        im = ax.imshow(dft.values, cmap=cmap, vmin=0.0, vmax=1.0, aspect=1)

        specie = species_name[layer]
        xlabel = '{specie:s} meiotic cell'.format(specie=specie)
        ax.set_xlabel(xlabel)
        ax.xaxis.set_label_position("top")

        specie = species_name[layer]
        ylabel = '{specie:s} soma'.format(specie=specie)
        ax.set_ylabel(ylabel)
        ax.yaxis.set_label_position("right")

        # Ticks
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        ax.set_yticklabels(yticklabels, rotation=0)
        ax.set_xticklabels(xticklabels, rotation=0)

    # Create colorbar
    axcb = fig.add_axes([0.050, 0.035, 0.12, 0.01])
    cbar = plt.colorbar(im, cax=axcb, orientation='horizontal', ticks=[0, 0.5, 1], boundaries=np.linspace(0, 1, 30))
    cbar.ax.tick_params(labelsize='small')
    cbar.ax.set_title('Similarity', rotation=0, fontsize='medium')



    plt.subplots_adjust(left=0.20, right=0.96, bottom=0.09, top=0.95, wspace=0.25, hspace=0.25)
    #plt.tight_layout()
    wIMGfile = 'images/module-proximity/img-module-proximity-{network:s}-{threshold:s}-combined.pdf'.format(celltype=celltype, network=network, threshold=threshold_str)
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()
