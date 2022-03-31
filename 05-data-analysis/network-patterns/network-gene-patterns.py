# coding=utf-8
# Author: Rion B Correia
# Date: July 22, 2020
#
# Description: Reads a threshold and a conserved multiLayer network to identify patterns .
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer, get_network_largest_connected_component
from itertools import combinations
from scipy.stats import ks_2samp
import matplotlib as mpl
from brokenaxes import brokenaxes
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt


def label_diff_x(ax, x0, x1, y, cap_size, y_text, text):
    xmid = x0 + ((x1 - x0) / 2)
    ax.plot([x0, x0, x1, x1], [y - cap_size, y, y, y - cap_size], color='black', lw=1)
    #
    ax.annotate(text, xy=(xmid, y), xytext=(xmid, y_text), zorder=10, ha='center', va='center', fontsize='medium')


def label_diff_y(ax, y0, y1, x, cap_size, x_text, text):
    ymid = y0 + ((y1 - y0) / 2)
    ax.plot([x + cap_size, x, x, x + cap_size], [y0, y0, y1, y1], color='black', lw=1)
    #
    ax.annotate(text, xy=(x, ymid), xytext=(x_text, ymid), zorder=10, ha='center', va='center', fontsize='medium')


def as_si(x, ndp=2):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


if __name__ == '__main__':

    celltypes = ['spermatocyte', 'enterocyte']
    #
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    layers = ['HS', 'MM', 'DM']

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
            'edgecolor': '#ffbb78'
        }
    }

    for celltype in celltypes:
        print("Loading celltype: {celltype:s}".format(celltype=celltype))
        #
        rGfile = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        G = nx.read_gpickle(rGfile)
        #
        for layer in layers:
            print('Separate layer {layer:s}'.format(layer=layer))
            Gl = get_network_layer(G, layer)
            #
            data[layer][celltype] = {
                'degree-centrality': [v for k, v in nx.degree_centrality(Gl).items()],
                'page-rank': [v for k, v in nx.pagerank(Gl, weight='weight').items()],
            }

    # Plot
    measures = ['degree-centrality', 'page-rank']
    titles = ['Degree Centrality', 'Page Rank']
    lylims = [((-0.5e-2, 3.1e-2), (3.9e-2, 4.7e-2)), ((-0.65e-4, 4.8e-4), (8.0e-4, 9.4e-4))]
    #
    #facecolors = ['#2ca02c', '#7f7f7f', '#ff7f0e']
    #edgecolors = ['#98df8a', '#c7c7c7', '#ffbb78']
    c_hs_sp = '#2ca02c'
    c_hs_en = '#bfe2bf'
    c_mm_sp = '#7f7f7f'
    c_mm_en = '#d8d8d8'
    c_dm_sp = '#ff7f0e'
    c_dm_en = '#ffd8b6'
    #
    flierprops = dict(marker='o', markersize=2, markeredgecolor='black', markeredgewidth=.05, rasterized=False)
    medianprops = dict(color='black')
    labels = ['HS Germ', 'HS Soma', 'MM Germ', 'MM Soma', 'DM Germ', 'DM Soma']

    for measure, ylims, title in zip(measures, lylims, titles):
        print('Plotting {measure:s}'.format(measure=measure))
        # Plot
        fig = plt.figure(figsize=(4.8, 3.6))
        bax = brokenaxes(ylims=ylims, hspace=0.05, despine=False)

        bax.set_title(title)
        bax.set_ylabel('Score', labelpad=25)
        """
        xvals = [0, 1, 2, 3, 4, 5]
        yvals = [
            data['HS']['spermatocyte'][measure],
            data['HS']['enterocyte'][measure],
            data['MM']['spermatocyte'][measure],
            data['MM']['enterocyte'][measure],
            data['DM']['spermatocyte'][measure],
            data['DM']['enterocyte'][measure],
        ]
        yvals_mean = [np.mean(vals) for vals in yvals]
        """
        xvals = [0, 1, 2.5, 3.5, 5.0, 6.0]
        yvals = [
            data['HS']['spermatocyte'][measure],
            data['HS']['enterocyte'][measure],
            data['MM']['spermatocyte'][measure],
            data['MM']['enterocyte'][measure],
            data['DM']['spermatocyte'][measure],
            data['DM']['enterocyte'][measure],
        ]
        
        xvals_hs = [xvals[0], xvals[1]]
        yvals_hs = [yvals[0], yvals[1]]
        yvals_mean_hs = [np.mean(vals) for vals in yvals_hs]
        #
        xvals_mm = [xvals[2], xvals[3]]
        yvals_mm = [yvals[2], yvals[3]]
        yvals_mean_mm = [np.mean(vals) for vals in yvals_mm]
        #
        xvals_dm = [xvals[4], xvals[5]]
        yvals_dm = [yvals[4], yvals[5]]
        yvals_mean_dm = [np.mean(vals) for vals in yvals_dm]
        #
        bp_hs_l = bax.boxplot(yvals_hs, positions=xvals_hs, labels=None, notch=True, patch_artist=True, widths=0.60,
            flierprops=flierprops, boxprops=dict(facecolor='#2ca02c'), medianprops=medianprops, showfliers=False, zorder=6)
        bax.plot(xvals_hs, yvals_mean_hs, marker='o', color='black', markersize=6, lw=0, zorder=8)

        bp_mm_l = bax.boxplot(yvals_mm, positions=xvals_mm, labels=None, notch=True, patch_artist=True, widths=0.60,
            flierprops=flierprops, boxprops=dict(facecolor='#2ca02c'), medianprops=medianprops, showfliers=False, zorder=6)
        bax.plot(xvals_mm, yvals_mean_mm, marker='o', color='black', markersize=6, lw=0, zorder=8)

        bp_dm_l = bax.boxplot(yvals_dm, positions=xvals_dm, labels=None, notch=True, patch_artist=True, widths=0.60,
            flierprops=flierprops, boxprops=dict(facecolor='#2ca02c'), medianprops=medianprops, showfliers=False, zorder=6)
        bax.plot(xvals_dm, yvals_mean_dm, marker='o', color='black', markersize=6, lw=0, zorder=8)
        #
        caps_l = [
            bp_hs_l[0]['caps'][0:2],
            bp_hs_l[0]['caps'][2:4],
            #
            bp_mm_l[0]['caps'][0:2],
            bp_mm_l[0]['caps'][2:4],
            #
            bp_dm_l[0]['caps'][0:2],
            bp_dm_l[0]['caps'][2:4],
        ]
        #
        for bp_hs in bp_hs_l:
            for patch, facecolor in zip(bp_hs['boxes'], [c_hs_sp, c_hs_en]):
                patch.set_facecolor(facecolor)
        for bp_mm in bp_mm_l:
            for patch, facecolor in zip(bp_mm['boxes'], [c_mm_sp, c_mm_en]):
                patch.set_facecolor(facecolor)
        for bp_dm in bp_dm_l:
            for patch, facecolor in zip(bp_dm['boxes'], [c_dm_sp, c_dm_en]):
                patch.set_facecolor(facecolor)

        # Left Tests
        y_hs_min = min(min(cap.get_ydata()) for cap in bp_hs_l[0]['caps'])
        y_mm_min = min(min(cap.get_ydata()) for cap in bp_mm_l[0]['caps'])
        Y_dm_min = min(min(cap.get_ydata()) for cap in bp_dm_l[0]['caps'])
        #
        y_hs_max = max(max(cap.get_ydata()) for cap in bp_hs_l[0]['caps'])
        y_mm_max = max(max(cap.get_ydata()) for cap in bp_mm_l[0]['caps'])
        y_dm_max = max(max(cap.get_ydata()) for cap in bp_dm_l[0]['caps'])

        ylimrange = (ylims[1][1] - ylims[0][0]) - (ylims[1][0] - ylims[0][1])
        text_padding = 0.04 * ylimrange
        #
        cap_size = 0.01 * ylimrange * -1
        for (idx_i, idx_j) in [(0, 1), (2, 3), (4, 5)]:
            x_i, x_j = xvals[idx_i], xvals[idx_j]
            ys_i, ys_j = yvals[idx_i], yvals[idx_j]
            caps = [cap for caps in caps_l[idx_i:idx_j + 1] for cap in caps]
            y = min([min(cap.get_ydata()) for cap in caps])
            y = y - (0.03 * ylimrange)
            y_text = y - text_padding
            stat, pvalue = ks_2samp(ys_i, ys_j, alternative='two-sided', mode='asymp')
            label_diff_x(ax=bax, x0=x_i, x1=x_j, y=y, cap_size=cap_size, y_text=y_text, text=r'${pvalue:s}$'.format(pvalue=as_si(pvalue)))

        cap_size = cap_size * -1.0
        y_last = None
        for (idx_i, idx_j) in [(0, 2), (2, 4), (0, 4)]:
            x_i, x_j = xvals[idx_i], xvals[idx_j]
            ys_i, ys_j = yvals[idx_i], yvals[idx_j]
            # first x position
            if y_last is None:
                caps = [cap for caps in caps_l[idx_i:idx_j + 1] for cap in caps]
                y = max([max(cap.get_ydata()) for cap in caps])
                y = y + (0.03 * ylimrange)
                y_last = y
            # next x position
            else:
                y = y_last
            y_text = y + (0.8 * text_padding)
            stat, pvalue = ks_2samp(ys_i, ys_j, alternative='two-sided', mode='asymp')
            label_diff_x(ax=bax, x0=x_i, x1=x_j, y=y, cap_size=cap_size, y_text=y_text, text=r'${pvalue:s}$'.format(pvalue=as_si(pvalue)))
            y_last += 0.08 * ylimrange

        # Legend
        pshs = mpl.patches.Patch(facecolor=c_hs_sp, edgecolor='k')
        psmm = mpl.patches.Patch(facecolor=c_mm_sp, edgecolor='k')
        psdm = mpl.patches.Patch(facecolor=c_dm_sp, edgecolor='k')
        #
        ls = plt.legend(handles=[pshs, psmm, psdm], labels=['', '', 'Spermatocyte'],
                        loc='upper left', ncol=3, handletextpad=0.5, handlelength=1.0, columnspacing=0.0)

        pehs = mpl.patches.Patch(facecolor=c_hs_en, edgecolor='k')
        pemm = mpl.patches.Patch(facecolor=c_mm_en, edgecolor='k')
        pedm = mpl.patches.Patch(facecolor=c_dm_en, edgecolor='k')
        le = plt.legend(handles=[pehs, pemm, pedm], labels=['', '', 'Enterocyte'],
                        loc='upper right', ncol=3, handletextpad=0.5, handlelength=1.0, columnspacing=0.0)

        bax.big_ax.add_artist(ls)
        bax.big_ax.add_artist(le)

        bax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)

        #bax.grid(axis='y')

        # Make big_ax the same limits
        bax.big_ax.set_xlim(bax.axs[0].get_xlim())
        ylim_min, ylim_max = bax.axs[0].get_ylim()
        bax.big_ax.set_ylim(ylim_min, ylim_max - ylim_min)

        # Disable duplicated "x10^-n" on axis
        bax.axs[1].yaxis.offsetText._visible = False
        #
        bax.axs[0].set_xticks([0.5, 3, 5.5])
        bax.axs[1].set_xticks([0.5, 3, 5.5])
        bax.axs[0].set_xticklabels(['Human', 'Mouse', 'Insect'])
        bax.axs[1].set_xticklabels(['Human', 'Mouse', 'Insect'])


        wIMGfile = 'images/network-{network:s}-{measure:s}.pdf'.format(network=network, measure=measure)
        ensurePathExists(wIMGfile)
        #
        #plt.subplots_adjust(left=0.05, right=0.98, bottom=0.10, top=0.90, wspace=0.4, hspace=0.25)
        #plt.tight_layout()
        plt.savefig(wIMGfile)
        plt.close()

