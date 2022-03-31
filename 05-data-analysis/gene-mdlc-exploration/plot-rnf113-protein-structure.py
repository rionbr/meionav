# coding=utf-8
# Author: Rion B Correia
# Date: Nov 24, 2021
#
# Description: 
#
#
import matplotlib as mpl
#
mpl.rc('font', size=10)  # controls default text sizes
mpl.rc('axes', titlesize=12)  # fontsize of the axes title
mpl.rc('axes', labelsize=12)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=10)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=10)  # fontsize of the tick labels
mpl.rc('legend', fontsize=10)  # legend fontsize
mpl.rc('figure', titlesize=12)  # fontsize of the figure title
#
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#

"""
RNF113B
Size: 1-322
Domains:
C3H1: 190-218; RING: 256-294

Proband’s (patient) RNF113B
Size: 1-304
Mutation: Position 186 (10nt deletion)
Domains:
None (due to frameshift)

dRNF113
Size: 1-357
Domains:
C3H1: 194-222; RING: 264-302
"""

if __name__ == '__main__':

    names = ['RNF113B', 'RNF113B del.', 'dRNF113']
    sizes = [322, 304, 357]
    listdomains = [
        [
            {'name': 'C3H1', 'loc': (190, 218), 'facecolor': '#ff9896', 'edgecolor': '#d62728'},
            {'name': 'RING', 'loc': (256, 294), 'facecolor': '#c5b0d5', 'edgecolor': '#9467bd'},
        ],
        [],
        [
            {'name': 'C3H1', 'loc': (194, 222), 'facecolor': '#ff9896', 'edgecolor': '#d62728'},
            {'name': 'RING', 'loc': (264, 302), 'facecolor': '#c5b0d5', 'edgecolor': '#9467bd'},
        ]
    ]

    # Plot
    fig = plt.figure(figsize=(3.3, 4.5))
    gs = gridspec.GridSpec(nrows=max(sizes), ncols=3)
    axes = []

    adjustable = 'datalim'
    aspect = 'auto'

    i = 0
    for name, size, domains in zip(names, sizes, listdomains):

        ax = plt.subplot(gs[:size, i])
        axes.append(ax)
        ax.set(adjustable=adjustable, aspect=aspect, anchor='NE')

        # Hide the right and top spines
        ax.xaxis.set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.set_title(name)
        ax.set(adjustable='box', aspect=1/35, anchor='NE')

        # Red background
        ax.fill_between([0.30, 0.70], 1, size, facecolor='#1f77b4', edgecolor='#353535', zorder=2)

        yticks = [1, 100]
        for domain in domains:
            name = domain.get('name')
            y0, y1 = domain.get('loc')
            cx, cy = 0.5, (y0 + (y1 - y0) / 2)
            facecolor = domain.get('facecolor')
            edgecolor = domain.get('edgecolor')
            #
            ax.fill_between([0, 1], y0, y1, facecolor=facecolor, edgecolor=edgecolor, zorder=3)
            ax.annotate(name, xy=(cx, cy), xytext=(cx, cy), xycoords='data', zorder=6, ha='center', va='center')
            #
            yticks.extend([y0, y1])

        if not len(domains):
            yticks += [200]
        yticks += [size]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)
        #
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-.1, size+1)
        ax.grid()
        ax.invert_yaxis()
        i += 1

    axes[1].fill_between([0.30, 0.70], 186, 304, facecolor='#aec7e8', edgecolor='#353535', zorder=3)
    axes[1].fill_between([0.30, 0.70], 186, 190, facecolor='#d62728', edgecolor='#353535', zorder=4)
    cx, cy = 1.02, 188
    axes[1].annotate('c.556_565del', xy=(cx, cy), xytext=(cx, cy), xycoords='data', zorder=6, ha='center', va='center', rotation=90)

    plt.subplots_adjust(left=0.04, right=0.93, bottom=0.04, top=0.94, wspace=0.3)

    wIMGFile = 'images/img-RNF113-protein-structure.pdf'
    plt.savefig(wIMGFile, dpi=300, pad_inches=0.0)
    plt.close()
