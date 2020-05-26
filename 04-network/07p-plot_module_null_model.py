# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads the similarity-celltype and plots the results.
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from utils import ensurePathExists
from data import *


# Confidence Interval
def calc_ci(x):
    mean, count, std = x.mean(), x.count(), x.std()
    return pd.Series({
        'mean': mean,
        'std': std,
        'ci95-max': (mean + 1.96 * std / math.sqrt(count)),
        'ci95-min': (mean - 1.96 * std / math.sqrt(count))})


def plot_module_null_model(celltype='spermatocyte', network='thr', threshold=0.5, layer='DM'):

    threshold_str = str(threshold).replace('.', 'p')
    print('Plotting {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))

    rCSVFile = 'results/module_null/{celltype:s}/module-null-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    df = pd.read_csv(rCSVFile, index_col=0, encoding='utf-8')

    df['mod-id'] = df['mod-id'].astype(str)
    df['values'] = df['values'].apply(eval)

    #colors = {'HS': '#2ca02c', 'MM': '#7f7f7f', 'DM': '#ff7f0e'}
    #color = colors[layer]

    data_cell = data_cells[celltype]
    modules = data_cell['modules-svd']['modules'][layer]


    fig = plt.figure(figsize=(6, 4))
    gs = gridspec.GridSpec(ncols=1, nrows=1, figure=fig)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_title('Gene (page)rank and SVD modules\n{celltype:s}-{network:s}-{threshold:.1f}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold, layer=layer))
    bins = np.linspace(0, 0.0005, 8)

    for mid, dfg in df.groupby('mod-id'):

        mname = dfg['mod-name'].iloc[0]
        print("M{mid:s} - {mname:s}".format(mid=str(mid), mname=mname))

        try:
            midcolor = int(mid)
        except:
            midcolor = mid
        else:
            midcolor = int(mid)
        for module in modules:
            if module['id'] == midcolor:
                color = module['facecolor']

        
        dfT = df.loc[((df['run'] == 'real') & (df['mod-id'] == mid)), :]

        dfN = df.loc[(~(df['run'] == 'real') & (df['mod-id'] == mid)), :]
        rN = []
        for idx, row in dfN.iterrows():
            dfNtmp = pd.DataFrame({'values': pd.Series(row['values']), 'run': int(row['run'])})
            rN.append(dfNtmp)
        dfNci = pd.concat(rN)
        dfNci = dfNci.sort_values(['run','values'], ascending=[True, False])
        dfNci['sequence'] = dfNci.groupby('run').cumcount()   
        dfNci = dfNci.pivot(index='sequence', columns='run', values='values')
        dfNci = dfNci.apply(calc_ci, axis='columns')


        real_vals = dfT['values'].tolist()[0]
        real_vals = sorted(real_vals, reverse=True)

        print(len(dfNci))
        print(len(real_vals))
        x = list(range(len(real_vals)))
        #hist_weights_real = np.ones_like(real_vals) / len(real_vals)
        #ax.hist(real_vals, bins=bins, density=False, weights=hist_weights_real, label="M{mid:s} - {mname:s}".format(mid=str(mid), mname=mname), edgecolor=color, facecolor=(0, 0, 0, 0), lw=1, zorder=5)
        ax.plot(x, real_vals, label="M{mid:s} - {mname:s}".format(mid=str(mid), mname=mname), color=color, marker='.', markersize=4, lw=2, zorder=5)
        
        
        #ax.plot(x, dfNci['mean'], color=color, marker=None, lw=1, zorder=5)
        ax.fill_between(x, dfNci['ci95-min'], dfNci['ci95-max'], color='black', alpha=0.4, zorder=3)
        # Ticks
        #xyticks = np.arange(len(index))
        #ax.set_xticks(xyticks)
        #ax.set_yticks(xyticks)

        # TickLabels
        #xyticklabels = ["M{cid:s}-{name:s}".format(cid=str(cid), name=name) for cid, name in index]
        #xyticklabels_fake = [x.split('-')[0] for x in xyticklabels]

        #ax.set_yticklabels(xyticklabels, rotation=0, fontsize='small')
        #ax.set_xticklabels(xyticklabels_fake, rotation=90, fontsize='small')


        # Create colorbar
        #cbar = plt.colorbar(im, cax=axcb, orientation='horizontal', ticks=[0, 0.5, 1], boundaries=np.linspace(0, 1, 30))
        #cbar.ax.tick_params(labelsize='small')
        #cbar.ax.set_title('Similarity', rotation=0, fontsize='small')

        #break

    # Null Model
    null_list = dfN['values'].tolist()
    null_vals = [i for l in null_list for i in l]
    hist_weights_null = np.ones_like(null_vals) / len(null_vals)
    #ax.hist(null_vals, bins=bins, density=False, weights=hist_weights_null, label='Null Model', edgecolor='gray', facecolor=(0, 0, 0, .10), lw=2, zorder=2)

    #ax.set_xticks(bins[0:-1:3])
    #ax.set_yticklabels(bins[0:-1:2])

    ax.set_ylabel('Gene (page)rank value')
    ax.set_xlabel('Rank')

    # Legend
    ax.legend(fontsize='small')

    ax.set_yscale('log')
    ax.grid()

    plt.subplots_adjust(left=0.10, right=0.96, bottom=0.13, top=0.86, wspace=0.0, hspace=0.0)
    wIMGfile = 'images/module_null/{celltype:s}/img-module-null-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer, mid=str(mid))
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()


if __name__ == '__main__':

    celltype = 'spermatocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    layer = 'HS'

    plot_module_null_model(celltype=celltype, network=network, threshold=threshold, layer=layer)

    asd
    for celltype in ['spermatocyte', 'enterocyte']:
        for layer in ['HS', 'MM', 'DM']:
            plot_module_null_model(celltype=celltype, network=network, threshold=threshold, layer=layer)

