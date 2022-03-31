# coding=utf-8
# Author: Rion B Correia
# Date: July 22, 2020
#
# Description: Reads a threshold and a conserved multiLayer network to identify network patterns.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
from itertools import combinations
from scipy.stats import ks_2samp
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
#
mpl.rc('font', size=10)  # controls default text sizes
mpl.rc('axes', titlesize=16)  # fontsize of the axes title
mpl.rc('axes', labelsize=12)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=10)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=10)  # fontsize of the tick labels
mpl.rc('legend', fontsize=10)  # legend fontsize
mpl.rc('figure', titlesize=16)  # fontsize of the figure title
#
mpl.rc('figure.subplot', left=0.15, right=0.9, bottom=0.11, top=0.88)
import matplotlib.pyplot as plt


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
    ldf = []
    for layer in layers:
        print("> Computing layer: {layer:s}".format(layer=layer))
        Gl = get_network_layer(Gt, layer=layer)

        dft = pd.DataFrame(index=Gl.nodes())

        # Layer
        dft['layer'] = layer
        # Conserved nodes
        dft['conserved'] = dft.index.map(nx.get_node_attributes(Gl, name='conserved'))
        dft['conserved'].fillna(False, inplace=True)
        # Degree Centrality
        dft['degree-centrality'] = dft.index.map(nx.degree_centrality(Gl))
        # Page Rank
        dft['page-rank'] = dft.index.map(nx.pagerank(Gl, weight='weight'))
        # append
        ldf.append(dft)

    df = pd.concat(ldf)

    # Separate Conserved / NonConserved
    dfc = df.loc[df['conserved'] == True, :].copy()
    dfn = df.loc[df['conserved'] == False, :].copy()


    flierprops = dict(marker='o', markersize=2, markeredgecolor='black', markeredgewidth=.5, rasterized=False)
    medianprops = dict(color='black')
    #colors = [('#2ca02c', '#d62728'), ('#7f7f7f', '#d62728'), ('#ff7f0e', '#d62728')]
    colors = [('#d62728', '#1f77b4'), ('#d62728', '#1f77b4'), ('#d62728', '#1f77b4')]
    measures = ['degree-centrality', 'page-rank']
    titles = ['Degree centrality', 'Page rank']

    for measure, title in zip(measures, titles):
        print("Plotting {measure:s}".format(measure=measure))

        # Plot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.1, 3.7))

        sc_hs = dfc.loc[(dfc['layer'] == 'HS'), measure]
        sc_mm = dfc.loc[(dfc['layer'] == 'MM'), measure]
        sc_dm = dfc.loc[(dfc['layer'] == 'DM'), measure]
        #
        sn_hs = dfn.loc[(dfn['layer'] == 'HS'), measure]
        sn_mm = dfn.loc[(dfn['layer'] == 'MM'), measure]
        sn_dm = dfn.loc[(dfn['layer'] == 'DM'), measure]
        #
        hs_xvals = [1, 2]
        hs_yvals = [sc_hs, sn_hs]
        mm_xvals = [4, 5]
        mm_yvals = [sc_mm, sn_mm]
        dm_xvals = [7, 8]
        dm_yvals = [sc_dm, sn_dm]
        #
        hs_yvals_mean = [v.mean() for v in hs_yvals]
        mm_yvals_mean = [v.mean() for v in mm_yvals]
        dm_yvals_mean = [v.mean() for v in dm_yvals]
        #
        bp_hs = ax.boxplot(hs_yvals, positions=hs_xvals, labels=['Conserved', 'Not conserved'], notch=True, patch_artist=True, widths=0.7, flierprops=flierprops, boxprops=dict(facecolor='#2ca02c'), medianprops=medianprops, showfliers=False, zorder=6)
        bp_mm = ax.boxplot(mm_yvals, positions=mm_xvals, labels=['Conserved', 'Not conserved'], notch=True, patch_artist=True, widths=0.7, flierprops=flierprops, boxprops=dict(facecolor='#7f7f7f'), medianprops=medianprops, showfliers=False, zorder=6)
        bp_dm = ax.boxplot(dm_yvals, positions=dm_xvals, labels=['Conserved', 'Not conserved'], notch=True, patch_artist=True, widths=0.7, flierprops=flierprops, boxprops=dict(facecolor='#ff7f0e'), medianprops=medianprops, showfliers=False, zorder=6)
        #
        ax.plot(hs_xvals, hs_yvals_mean, marker='o', color='black', markersize=6, lw=0, zorder=8)
        ax.plot(mm_xvals, mm_yvals_mean, marker='o', color='black', markersize=6, lw=0, zorder=8)
        ax.plot(dm_xvals, dm_yvals_mean, marker='o', color='black', markersize=6, lw=0, zorder=8)
        #
        ax.set_title(title)
        ax.set_ylabel('Score')
        ax.set_xticks([1.5, 4.5, 7.5])
        ax.set_xticklabels(['Human', 'Mouse', 'Insect'])
        #
        #ax.grid(axis='y')
        #
        y_hs = max(max(cap.get_ydata()) for cap in bp_hs['caps'])
        y_mm = max(max(cap.get_ydata()) for cap in bp_mm['caps'])
        y_dm = max(max(cap.get_ydata()) for cap in bp_dm['caps'])
        y_species = max(y_hs, y_mm, y_dm)
        #
        y_hs_start = 0.05 * y_hs
        y_mm_start = 0.05 * y_mm
        y_dm_start = 0.05 * y_dm
        cap_size = 0.02 * y_species
        y_text = text_padding = 0.05 * y_species
        #
        stat, pvalue = ks_2samp(sc_hs, sn_hs, alternative='two-sided', mode='asymp')
        #  print("KS test: stat={stat:.3f}, p-value={pvalue:.2e}".format(stat=stat, pvalue=pvalue))
        y = y_hs + y_hs_start
        y_text = y + text_padding
        text = '****' if pvalue <= 0.001 else pvalue  # r'${pvalue:s}$'.format(pvalue=as_si(pvalue))
        label_diff(ax=ax, x0=1, x1=2, y=y, cap_size=cap_size, y_text=y_text, text=text)
        #
        stat, pvalue = ks_2samp(sc_mm, sn_mm, alternative='two-sided', mode='asymp')
        #  print("KS test: stat={stat:.3f}, p-value={pvalue:.2e}".format(stat=stat, pvalue=pvalue))
        y = y_mm + y_mm_start
        y_text = y + text_padding
        text = '****' if pvalue <= 0.001 else pvalue  # r'${pvalue:s}$'.format(pvalue=as_si(pvalue))
        label_diff(ax=ax, x0=4, x1=5, y=y, cap_size=cap_size, y_text=y_text, text=text)

        stat, pvalue = ks_2samp(sc_dm, sn_dm, alternative='two-sided', mode='asymp')
        #  print("KS test: stat={stat:.4f}, p-value={pvalue:.2e}".format(stat=stat, pvalue=pvalue))
        y = y_dm + y_dm_start
        y_text = y + text_padding
        text = '****' if pvalue <= 0.001 else pvalue  # r'${pvalue:s}$'.format(pvalue=as_si(pvalue))
        label_diff(ax=ax, x0=7, x1=8, y=y, cap_size=cap_size, y_text=y_text, text=text)

        for bp, (color_con, color_notcon) in zip([bp_hs, bp_mm, bp_dm], colors):
            for patch, color in zip(bp['boxes'], [color_con, color_notcon]):
                patch.set_facecolor(color)

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)

        xlim_min, xlim_max = ax.get_ylim()
        xlim_max = 1.04 * xlim_max
        ax.set_ylim((xlim_min, xlim_max))
        ax.set_xlim(0.0, 9)
        #
        psc = mpl.patches.Patch(facecolor='#d62728', edgecolor='k')
        psn = mpl.patches.Patch(facecolor='#1f77b4', edgecolor='k')
        #ax.legend((psc, psn), ('Conserved', 'Non-conserved'),
        #    loc='upper right',
        #    ncol=2, handletextpad=0.5, handlelength=1.25, columnspacing=1.25)
        #

        wIMGfile = 'images/network-{network:s}-{measure:s}.pdf'.format(network='thr-conserved', measure=measure)
        ensurePathExists(wIMGfile)
        #
        #plt.subplots_adjust(left=0.03, right=0.90, bottom=0.12, top=0.87, wspace=0.0, hspace=0.0)
        #plt.tight_layout()
        plt.savefig(wIMGfile)
        plt.close()

