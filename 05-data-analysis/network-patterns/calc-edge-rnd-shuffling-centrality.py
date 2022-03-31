# coding=utf-8
# Author: Rion B Correia
# Date: Jan 31, 2021
#
# Description: null model for - Reads a threshold and a conserved multiLayer network to identify network patterns.
#
#
import random
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
from itertools import combinations
from scipy.stats import ks_2samp, ttest_ind
import matplotlib as mpl
import matplotlib.ticker as mtick
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
#
mpl.rc('font', size=10)  # controls default text sizes
mpl.rc('axes', titlesize=14)  # fontsize of the axes title
mpl.rc('axes', labelsize=12)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=10)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=10)  # fontsize of the tick labels
mpl.rc('legend', fontsize=10)  # legend fontsize
mpl.rc('figure', titlesize=16)  # fontsize of the figure title
#
mpl.rc('figure.subplot', left=0.19, right=0.95, bottom=0.15, top=0.85)
import matplotlib.pyplot as plt


def generate_random_nonedge(G):
    while True:
        i, j = np.random.choice(G.nodes(), 2)
        if not G.has_edge(i, j):
            break
    return (i, j)

"""
def label_diff(ax, xs, xe, text, bp):
    caps = bp['caps']
    y = 1.05 * max([max(cap.get_ydata()) for cap in caps])
    ydip = 0.98 * y
    ytext = 1.04 * y

    xmid = xs + ((xe - xs) / 2)
    ax.plot([xs, xs, xe, xe], [ydip, y, y, ydip], color='black', lw=1)
    #lde = lines.Line2D([xs, ydip], [xs, y])
    #lmi = lines.Line2D([xs, y], [xe, y])
    #ldr = lines.Line2D([xe, y], [xe, ydip])
    #ax.lines.extend([lde, lmi, ldr])
    ax.annotate(text, xy=(xmid, y), xytext=(xmid, ytext), zorder=10, ha='center', va='center', fontsize='x-small')
"""

def label_diff(ax, x0, x1, y, cap_size, y_text, text):
    xmid = x0 + ((x1 - x0) / 2)
    ax.plot([x0, x0, x1, x1], [y-cap_size, y, y, y-cap_size], color='black', lw=1)
    #
    ax.annotate(text, xy=(xmid, y), xytext=(xmid, y_text), zorder=10, ha='center', va='center', fontsize='medium')


def as_si(x, ndp=2):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


if __name__ == '__main__':

    celltype = 'spermatocyte'

    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    layers = ['HS', 'MM', 'DM']
    layers_str = {
        'HS': 'Human',
        'MM': 'Mouse',
        'DM': 'Insect'
    }

    rGtfile = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='thr', threshold=threshold_str)
    Gt = nx.read_gpickle(rGtfile)

    rGcfile = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network='conserved', threshold=threshold_str)
    Gc = nx.read_gpickle(rGcfile)

    # Remove Cross edges
    Gt.remove_edges_from([(i, j) for i, j, d in Gt.edges(data=True) if d['type'] == 'cross'])

    # Identify conserved nodes
    nx.set_node_attributes(Gt, values={k: True for k in Gc.nodes()}, name='conserved')

    data = {
        'HS': {
            'name': 'Human',
            'facecolor': '#2ca02c',
            'edgecolor': '#98df8a',
        },
        'MM': {
            'name': 'Mouse',
            'facecolor': '#7f7f7f',
            'edgecolor': '#c7c7c7',
        },
        'DM': {
            'name': 'Insect',
            'facecolor': '#ff7f0e',
            'edgecolor': '#ffbb78',
        }
    }

    #
    measures = ['degree-centrality', 'page-rank']

    r = []
    for layer in layers:
        print("> Computing layer: {layer:s}".format(layer=layer))
        Gl = get_network_layer(Gt, layer=layer)
        #
        nodes = list(Gl.nodes())
        edges = list(Gl.edges())
        #
        dfn = pd.DataFrame(index=Gl.nodes())
        dfe = nx.to_pandas_edgelist(Gl)
        dfe = dfe.loc[:, ['source', 'target', 'weight']]
        #
        print("> Generating non-edges")
        dfne = pd.DataFrame(nx.non_edges(Gl), columns=['source', 'target'])

        # Conserved
        dfn['conserved'] = dfn.index.map(nx.get_node_attributes(Gl, name='conserved'))
        dfn['conserved'].fillna(False, inplace=True)

        # Iterate percentages of swapping
        for percentage in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            print("Percentage: {percentage:.2f}".format(percentage=percentage))
            print("> Copy G")
            Glt = Gl.copy()

            if percentage > 0.0:
                print("> Edge swaping")
                n_e_to_swap = int(percentage * Glt.number_of_edges())

                dfe_chosen = dfe.sample(n=n_e_to_swap, replace=False).copy()
                dfne_chosen = dfne.sample(n=n_e_to_swap, replace=False).copy()
                # remove
                Glt.remove_edges_from(dfe_chosen[['source', 'target']].values.tolist())

                dfe_chosen['source'] = dfne_chosen['source'].tolist()
                dfe_chosen['target'] = dfne_chosen['target'].tolist()

                # add
                records = dfe_chosen.set_index(['source', 'target']).to_dict(orient='index')
                tuples = [(i, j, d) for (i, j), d in records.items()]
                Glt.add_edges_from(tuples)

            print("> Calculating measures")
            # Degree Centrality
            dfn['degree-centrality'] = dfn.index.map(nx.degree_centrality(Glt))
            # Page Rank
            dfn['page-rank'] = dfn.index.map(nx.pagerank(Glt, weight='weight'))

            for measure in measures:

                sc = dfn.loc[(dfn['conserved'] == True), measure]
                sn = dfn.loc[(dfn['conserved'] == False), measure]

                #stat, pvalue = ks_2samp(sc, sn, alternative='two-sided', mode='auto')
                #stat, pvalue = ttest_ind(sc, sn, equal_var=False) #, alternative='two-sided', trim=0)
                #stat, pvalue = mannwhitneyu(sc, sn, alternative="two-sided")

                r.append((layer, percentage, measure, sc.mean(), sn.mean(), sc.std(), sn.std(), sc.median(), sn.median()))
    
    #columns = ['layer', 'percentage', 'measure', 'stat', 'p-value']
    columns = ['layer', 'percentage', 'measure', 'sc-mean', 'sn-mean', 'sc-std', 'sn-std', 'sc-median', 'sn-median']
    df = pd.DataFrame(r, columns=columns)

    # Export
    df.to_csv('results/csv-null-model.csv'.format(layer=layer))

    #colors = [('#2ca02c', '#d62728'), ('#7f7f7f', '#d62728'), ('#ff7f0e', '#d62728')]
    colors = [('#d62728', '#1f77b4'), ('#d62728', '#1f77b4'), ('#d62728', '#1f77b4')]
    cols = ['degree-centrality', 'page-rank']
    titles = ['Degree centrality', 'Page rank']

    for col, measure, title in zip(cols, measures, titles):
        print("Plotting {measure:s}".format(measure=measure))

        # Plot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.1, 3.7))

        scn_hs = df.loc[(df['layer'] == 'MM') & (df['measure'] == col), :]
        #scn_mm = df.loc[(df['layer'] == 'MM'), col]
        #scn_dm = df.loc[(df['layer'] == 'DM'), col]
        #
        pc_hs = ax.errorbar(x=scn_hs['percentage'], y=scn_hs['sc-mean'], yerr=scn_hs['sc-std'],
            color='#d62728', ecolor='#ff9896', elinewidth=2, marker='o', capsize=4, label='Conserved', zorder=6)
        pn_hs = ax.errorbar(x=scn_hs['percentage'], y=scn_hs['sn-mean'], yerr=scn_hs['sn-std'],
            color='#1f77b4', ecolor='#aec7e8', elinewidth=2, marker='o', capsize=4, label='Not Conserved', zorder=5)
    
        #
        layer_str = layers_str[layer]
        title = '{layer:s} network rewiring\n{title:s}'.format(layer=layer_str, title=title)
        ax.set_title(title)
        ax.set_ylabel('Mean (SD)')
        ax.set_xlabel('Edge rewiring')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
        #ax.set_xticks(scn_hs['percentage'])
        ax.xaxis.set_major_formatter(mtick.PercentFormatter(1.0))
        #
        ax.legend()
        #
        wIMGfile = 'images/network-{network:s}-{layer:s}-null-{measure:s}.pdf'.format(network='thr-conserved', layer=layer, measure=measure)
        ensurePathExists(wIMGfile)
        #
        #plt.subplots_adjust(left=0.03, right=0.90, bottom=0.12, top=0.87, wspace=0.0, hspace=0.0)
        #plt.tight_layout()
        plt.savefig(wIMGfile)
        plt.close()

