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


def plot_module_proximity(celltype='spermatocyte', network='thr', threshold=0.5):

    threshold_str = str(threshold).replace('.', 'p')

    #
    print('Loading .csv')
    rCSVFile = 'results/proximity/svd-proximity-{celltype:s}-{network:s}-{threshold:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str)

    df = pd.read_csv(rCSVFile, index_col=0, encoding='utf-8')
    print(df)

    if celltype == 'spermatocyte':
        df = df.loc[((df['id-i']<=9) & (df['id-j']<=9)), :]

    # Index
    dfi = df[['id-i', 'name-i']].drop_duplicates().rename(columns={'id-i': 'id', 'name-i': 'name'})
    dfj = df[['id-j', 'name-j']].drop_duplicates().rename(columns={'id-j': 'id', 'name-j': 'name'})
    dfij = pd.concat([dfi, dfj]).drop_duplicates().sort_values('id').reset_index(drop=True)
    index = dfij.set_index(['id', 'name']).index

    print('Plot')
    fig = plt.figure(figsize=(8, 2.6))
    gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)
    ax_hs_mm = fig.add_subplot(gs[0, 0])
    ax_hs_dm = fig.add_subplot(gs[0, 1])
    ax_mm_dm = fig.add_subplot(gs[0, 2])

    axcb = fig.add_axes([0.06, 0.12, 0.14, 0.021])

    title = 'Jaccard similarity between {celltype:s} SVD modules'.format(celltype=celltype)
    fig.suptitle(title, )

    #
    layer_pairs = [('HS', 'MM'), ('HS', 'DM'), ('MM', 'DM')]
    axes = [ax_hs_mm, ax_hs_dm, ax_mm_dm]
    species_name = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    cmap = cm.get_cmap('jet')
    cmap.set_bad(color='#bdbdbd')

    for (layer_i, layer_j), ax in zip(layer_pairs, axes):
        #
        dft = df.loc[(df['layer-i'] == layer_i) & (df['layer-j'] == layer_j), :]
        dft = dft.pivot_table(index=['id-i', 'name-i'], columns=['id-j', 'name-j'], values='proximity', aggfunc='first')
        dft = dft.reindex(index=index, columns=index, fill_value=np.nan)
        dft = dft.T
        #
        species_i = species_name[layer_i]
        species_j = species_name[layer_j]
        # ax.set_title('{species_i:s} vs {species_j:s}'.format(species_i=species_i, species_j=species_j))

        im = ax.imshow(dft.values, cmap=cmap, vmin=0.0, vmax=1.0)

        ax.set_xlabel(species_i)
        ax.xaxis.set_label_position("top")
        ax.set_ylabel(species_j)
        ax.yaxis.set_label_position("right")

    # Ticks
    xyticks = np.arange(len(index))
    ax_hs_mm.set_xticks(xyticks)
    ax_hs_mm.set_yticks(xyticks)
    #
    ax_hs_dm.set_xticks(xyticks)
    ax_hs_dm.set_yticks(xyticks)
    #
    ax_mm_dm.set_xticks(xyticks)
    ax_mm_dm.set_yticks(xyticks)

    # TickLabels
    xyticklabels = ["M{cid:s}-{name:s}".format(cid=str(cid), name=name) for cid, name in index]
    xyticklabels_fake = [x.split('-')[0] for x in xyticklabels]
    
    ax_hs_mm.set_yticklabels(xyticklabels, rotation=0, fontsize='small')
    ax_hs_mm.set_xticklabels(xyticklabels_fake, rotation=90, fontsize='small')
    #
    ax_hs_dm.set_yticklabels(xyticklabels_fake, rotation=0, fontsize='small')
    ax_hs_dm.set_xticklabels(xyticklabels_fake, rotation=90, fontsize='small')
    #
    ax_mm_dm.set_yticklabels(xyticklabels_fake, rotation=0, fontsize='small')
    ax_mm_dm.set_xticklabels(xyticklabels_fake, rotation=90, fontsize='small')

    # Create colorbar
    cbar = plt.colorbar(im, cax=axcb, orientation='horizontal', ticks=[0, 0.5, 1], boundaries=np.linspace(0, 1, 30))
    cbar.ax.tick_params(labelsize='small')
    cbar.ax.set_title('Similarity', rotation=0, fontsize='small')

    plt.subplots_adjust(left=0.25, right=0.96, bottom=0.08, top=0.98, wspace=0.50, hspace=0.0)
    wIMGfile = 'images/proximity/img-svd-proximity-{celltype:s}-{network:s}-{threshold:s}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str)
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()


if __name__ == '__main__':

    celltype = 'spermatocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5

    # plot_module_proximity(celltype=celltype, network=network, threshold=threshold)

    for celltype in ['spermatocyte', 'enterocyte']:
        plot_module_proximity(celltype=celltype, network=network, threshold=threshold)
