# coding=utf-8
# Author: Rion B Correia
# Date: Aug 06, 2019
#
# Description: Plots results of screened DM genes
#
# Instructions:
#
import numpy as np
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt


if __name__ == '__main__':

    fig = plt.figure(figsize=(4, 4))

    #
    # FPKM
    #
    norm_fpkm = mpl.colors.Normalize(vmin=0, vmax=11.79)
    cmap_fpkm = mpl.cm.Reds
    ax_fpkm_cb = fig.add_axes([0.5, 0.6, 0.32, 0.03])  # x, y, width, height
    ax_fpkm_cb.set_title('FPKM', fontsize='small', loc='left')
    cb_fpkm = mpl.colorbar.ColorbarBase(ax=ax_fpkm_cb, cmap=cmap_fpkm, norm=norm_fpkm, orientation='horizontal')
    cb_fpkm.set_label(r'log$_2$(FPKM+1)', fontsize='small')
    cb_fpkm.ax.tick_params(labelsize='small')

    #
    # RNAi
    #
    ax_rnai_cb = fig.add_axes([0.1, 0.6, 0.03, (0.05 * 2)])  # x, y, width, height
    ax_rnai_cb.set_title('Validated RNAi', fontsize='small', loc='left')
    cmap_rnai = mpl.colors.ListedColormap(['white', '#17becf'])
    bounds_rnai = [0, 1, 2]
    ticks_rnai = [0.5, 1.5]
    norm_rnai = mpl.colors.BoundaryNorm(bounds_rnai, cmap_rnai.N)
    cb_rnai = mpl.colorbar.ColorbarBase(ax=ax_rnai_cb, cmap=cmap_rnai, norm=norm_rnai, boundaries=bounds_rnai, extend='neither', ticks=ticks_rnai, spacing='uniform', orientation='vertical')
    cb_rnai.set_ticklabels(['No', 'Yes'])
    cb_rnai.ax.set_xticklabels(['No', 'Yes'], rotation=90, fontsize='small')
    cb_rnai.ax.tick_params(labelsize='small')

    #
    # Code
    #
    code_label = {
        'A': 'Meiotic',
        'B': 'Post-meiotic',
        'C': 'Gametes',
        'D': 'Pre-meiotic',
        'E': 'General impairment of spermatogenesis',
        'F': 'Undetectable',
        'G': 'Unspecified ',
        'H': 'Non-germ cell autonomous'
    }
    code_color = {
        'A': '#d62728',
        'B': '#ce6dbd',
        'C': '#756bb1',
        'D': '#c7e9c0',
        'E': '#9edae5',
        'F': '#fdd0a2',
        'G': '#dadaeb',
        'H': '#bdbdbd'
    }

    ax_code_cb = fig.add_axes([0.1, 0.05, 0.03, (0.05 * 8)])  # x, y, width, height
    ax_code_cb.set_title('Phenotype Code', fontsize='small', loc='left')
    cmap_code = mpl.colors.ListedColormap(list(code_color.values())[::-1])
    bounds_code = np.arange(len(code_label) + 1)
    ticks_code = np.arange(len(code_label) + 1) + 0.5
    norm_code = mpl.colors.BoundaryNorm(bounds_code, cmap_code.N)
    cb_code = mpl.colorbar.ColorbarBase(ax=ax_code_cb, cmap=cmap_code, norm=norm_code, boundaries=bounds_code, extend='neither', ticks=ticks_code, spacing='uniform', orientation='vertical')
    cb_code.set_ticklabels(list(code_label.values())[::-1])
    cb_code.ax.set_xticklabels(list(code_label.values()), rotation=90, fontsize='small')
    cb_code.ax.tick_params(labelsize='small')

    # Save Figure
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.05, top=0.95, wspace=0.2, hspace=0.2)
    file = 'images/img-core_DM-screened-legend.pdf'
    fig.savefig(file)
