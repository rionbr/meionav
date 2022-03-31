# coding=utf-8
# Author: Rion B Correia
# Date: Aug 06, 2019
#
# Description: Plots results of screened DM genes
#
# Instructions:
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
from matplotlib import colors
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from collections import OrderedDict


def value_to_color(x, cmap, norm):
    rgb = cmap(norm(x))[:3]
    return colors.rgb2hex(rgb)


def calc_control_mean_std_fert_rate(x):
    fertrate = x['hatched'] / x['eggs']
    return pd.Series({'mean fert-rate': fertrate.mean(), 'std fert-rate': fertrate.std()})


if __name__ == '__main__':

    # Load genes
    df = pd.read_csv('../02-core_genes/results/pipeline-core/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])

    # Load Screened data
    dfs = pd.read_csv('data/core_DM_screened_2021-06-07.csv', index_col=0)

    # Load Control data
    dfc = pd.read_csv('data/screened_DM_controls.csv', index_col=0)
    dfc = dfc.groupby(dfc.index).apply(calc_control_mean_std_fert_rate)

    # Load FPKM data
    dfFPKM = pd.read_csv('../02-core_genes/results/FPKM/DM/DM-FPKM-spermatocyte.csv.gz', index_col=0, usecols=['id_gene', 'FPKM'])

    # Load Module data
    dfM = pd.read_csv('../05-data-analysis/entropy-based-modules/results/pca-entropy/spermatocyte/DM/pca-spermatocyte-thr-0p5-DM-modules.csv.gz', index_col=0)
    # remove duplicates, keep first module
    dfM = dfM.loc[~dfM.index.duplicated(keep='first'), :]

    dfs_only = dfs.loc[~dfs['FT1 eggs'].isnull(), :]

    #status_cats = ['Screened', 'To be crossed', 'Pending', 'Reorder']
    #dfs['Status'] = pd.Categorical(dfs['Status'], categories=status_cats, ordered=True)
    #df['Status'] = dfs['Status']

    cols = ['FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched']
    df[cols] = dfs_only[cols]

    # Calculations
    df['total-eggs'] = 0
    df['total-hatched'] = 0
    for ft in range(1, 5):

        col_eggs = 'FT{:d} eggs'.format(ft)
        col_hatched = 'FT{:d} hatched'.format(ft)
        col_fertate = 'FT{:d} fert-rate'.format(ft)
        df[col_fertate] = df[col_hatched] / df[col_eggs]
        df['total-eggs'] += df[col_eggs]
        df['total-hatched'] += df[col_hatched]

    # Mean/SD
    df['mean fert-rate'] = df[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].mean(axis=1)
    df['std fert-rate'] = df[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].std(axis=1)

    # print(dfs.head())
    # print(dfs.loc[dfs['MM pheno code'].str.len() > 1, :])
    df['RNAi'] = dfs['Previous ref to RNAi working?']
    df['our-DM-code'] = dfs['Our DM pheno code']
    df['ext-DM-code'] = dfs['Others DM pheno code']
    df['ext-MM-code'] = dfs['MM pheno code']
    df['ext-HS-code'] = dfs['HS pheno code']

    #Module
    df['module-id'] = dfM['module-id']
    # Print/Export only matches with modules
    df.loc[df['module-id'].notnull(), ['gene', 'module-id']].sort_values('module-id').to_csv('results/csv-screened-genes-in-modules.csv')

    # Only those with a phenotype
    df = df.loc[df['our-DM-code'].notnull(), :]

    # FPKM
    df['FPKM'] = dfFPKM['FPKM']
    df['logFPKM'] = df['FPKM'].apply(lambda x: np.log2(x + 1))

    df = df.sort_values(['mean fert-rate', 'module-id', 'gene'], ascending=[True, True, True]).reset_index()
    df.index += 1

    maxfpkm, minfpkm = df['logFPKM'].max(), df['logFPKM'].min()

    code_label = OrderedDict({
        'D': 'Pre-meiotic',
        'A': 'Meiotic',
        'B': 'Post-meiotic',
        'C': 'Gametes',
        'F': 'Undetectable',
        # space
        'E': 'General',
        'H': 'Non-germ cell',
        'G': 'Unspecified ',
    })

    code_color = OrderedDict({
        'D': 'green',
        'A': 'red',
        'B': 'mediumblue',
        'C': 'darkmagenta',
        'F': 'gold',
        # space
        'E': 'cyan',
        'H': 'dimgray',
        'G': 'lightgray',
    })

    module_color = OrderedDict({
        1: '#e377c2',
        2: '#2ca02c',
        3: '#ff7f0e',
        4: '#8c564b',
        5: '#9467bd',
        6: '#bcbd22',
        11: '#6b6ecf',
        12: '#8ca252',
        np.nan: '#1f77b4',
    })
    module_label = OrderedDict({
        1: 'M1-Ubiquitination',
        2: 'M2-Splicing',
        3: 'M3-Translation',
        4: 'M4-rRNA regulation',
        5: 'M5-Vesicle transport',
        6: 'M6-Respiration',
        11: 'M11-Metabolism',
        12: 'M12-Pep. dephospho.',
        np.nan: 'Unspecified'
    })

    # print(df.head())
    # print(df.tail())
    print("logFPKM: {:.2f}/{:.2f}".format(minfpkm, maxfpkm))

    #
    # Plot Page 1
    #
    print("Plotting")
    n_per_page = 100
    dfg = df.groupby(np.arange(len(df)) // n_per_page)
    number_of_pages = len(dfg)
    for page, dft in dfg:
        # Debug
        #if page < 2:
        #    continue
        #
        n_this_page = len(dft)
        if (page + 1) == number_of_pages:
            last_page = True
            n_this_page += 5  # add the 3 + 2 = 5 (3 for manual fix; 2 controls)
        else:
            last_page = False
        print("> Page: {:d} of {:d}".format(page + 1, number_of_pages))
        print("> Points in this page: {:d}".format(n_this_page))

        fig = plt.figure(figsize=(6, 3.0))
        # fig.suptitle('Core metazoan meiotic genes'.format(page, number_of_pages))

        gs = gridspec.GridSpec(nrows=16, ncols=n_per_page)
        if not last_page:
            n_for_grid_width = n_this_page
        else:
            n_for_grid_width = n_this_page - 4
        # ax_rnai = plt.subplot(gs[1, :n_for_grid_width])
        ax_fpkm = plt.subplot(gs[1:2, :n_for_grid_width])
        ax_ext_hs = plt.subplot(gs[2:3, :n_for_grid_width])
        ax_ext_mm = plt.subplot(gs[3:4, :n_for_grid_width])
        ax_ext_dm = plt.subplot(gs[4:5, :n_for_grid_width])
        ax_our_dm = plt.subplot(gs[6:7, :n_for_grid_width])
        ax_fert = plt.subplot(gs[7:15, :n_for_grid_width])
        if last_page:
            ax_ctr = plt.subplot(gs[7:15, n_for_grid_width + 2:n_for_grid_width + 7])

        adjustable = 'datalim'
        aspect = 'auto'
        #ax_rnai.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_fpkm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_ext_hs.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_ext_mm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_ext_dm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_our_dm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_fert.set(adjustable=adjustable, aspect=aspect, anchor='NE')

        norm_fpkm = mpl.colors.Normalize(vmin=minfpkm, vmax=maxfpkm)
        cmap_fpkm = mpl.cm.Reds

        rotation = 50
        s = 7
        marker = 's'
        n = len(df)
        xticks = list(np.arange(dft.index.min() - 1, dft.index.max() + 1, 20))
        xticklabels = xticks
        xlimmin, xlimmax = dft.index.min() - 1, dft.index.max() + 1
        #
        # RNAi
        #
        """
        data_rnai = dft.loc[dft['RNAi'] == 'Yes', 'RNAi']
        x = data_rnai.index
        y = np.zeros(len(x))
        c = '#17becf'

        sc_rnai = ax_rnai.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_rnai.set_yticks([0])
        ax_rnai.set_yticklabels([''])
        ax_rnai.set_xticklabels([])
        ax_rnai.set_ylim(-0.2, 0.2)
        ax_rnai.set_xlim(-1, len(dft))
        """

        #
        # Expression (FPKM)
        #
        x = dft['logFPKM'].index
        y = np.zeros(len(x))
        c = dft['logFPKM'].apply(value_to_color, args=(cmap_fpkm, norm_fpkm))
        sc_fpkm = ax_fpkm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_fpkm.set_yticks([0])
        ax_fpkm.set_yticklabels(['FPKM'])
        ax_fpkm.set_xticks(xticks)
        ax_fpkm.set_xticklabels([])
        ax_fpkm.set_ylim(-0.5, 0.5)  # Adjusting this makes the plot shrink
        ax_fpkm.set_xlim(xlimmin, xlimmax)
        #ax_fpkm.set_xlim(-1, len(dft) + 1)

        #
        # External HS Phenotype
        #
        data_ext_hs = dft.loc[~dft['ext-HS-code'].isnull(), 'ext-HS-code']
        x = data_ext_hs.index
        y = np.zeros(len(x))
        c = data_ext_hs.map(code_color)

        sc_ext_hs = ax_ext_hs.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_ext_hs.set_yticks([0])
        ax_ext_hs.set_yticklabels(['HS data'])
        ax_ext_hs.set_xticklabels([])
        ax_ext_hs.set_ylim(-0.5, 0.5)
        ax_ext_hs.set_xlim(xlimmin, xlimmax)

        #
        # External MM Phenotype
        #
        # (these lines solve the problem when an MM phenotype has two codes, e.g., A/B)
        data_ext_mm = dft.loc[~dft['ext-MM-code'].isnull(), 'ext-MM-code']
        if len(data_ext_mm):
            data_tmp = data_ext_mm.str.split('/').apply(pd.Series)
            data_tmp = pd.melt(data_tmp.reset_index(), id_vars='index', value_vars=data_tmp.columns.tolist())
            data_tmp = data_tmp.set_index('index').dropna(subset=['value'])
            data_tmp.loc[(data_tmp.index.duplicated(keep=False) & (data_tmp['variable'] == 0)), 'variable'] = -0.2
            data_tmp.loc[(data_tmp.index.duplicated(keep=False) & (data_tmp['variable'] == 1)), 'variable'] = +0.2
            #
            x = data_tmp.index.values
            y = data_tmp['variable'].values
            c = data_tmp['value'].map(code_color).values
            sc_ext_mm = ax_ext_mm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_ext_mm.set_yticks([0])
        ax_ext_mm.set_yticklabels(['MM data'])
        ax_ext_mm.set_xticklabels([])
        ax_ext_mm.set_ylim(-0.5, 0.5)
        ax_ext_mm.set_xlim(xlimmin, xlimmax)

        #
        # External DM Phenotype
        #
        data_ext_dm = dft.loc[~dft['ext-DM-code'].isnull(), 'ext-DM-code']
        x = data_ext_dm.index
        y = np.zeros(len(x))
        c = data_ext_dm.map(code_color)

        sc_ext_dm = ax_ext_dm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_ext_dm.set_yticks([0])
        ax_ext_dm.set_yticklabels(['DM data'])

        ax_ext_dm.set_xticklabels([])
        ax_ext_dm.set_ylim(-0.5, 0.5)
        ax_ext_dm.set_xlim(xlimmin, xlimmax)
        # ax_ext_dm.grid(axis='y', linewidth=0.5)

        #
        # Our DM Phenotype
        #
        data_our_dm = dft.loc[~dft['our-DM-code'].isnull(), 'our-DM-code']
        x = data_our_dm.index
        y = np.zeros(len(x))
        c = data_our_dm.map(code_color)

        sc_our_dm = ax_our_dm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_our_dm.set_yticks([0])
        ax_our_dm.set_yticklabels(['DM pheno.'])
        ax_our_dm.set_xticklabels([])
        ax_our_dm.set_ylim(-0.2, 0.2)
        ax_our_dm.set_xlim(xlimmin, xlimmax)
        # ax_our_dm.grid(axis='y', linewidth=0.5)

        #
        # DM Expression Values
        #
        xs = dft.index
        ys = dft['mean fert-rate']
        es = dft['std fert-rate']
        cs = dft['module-id'].map(module_color)
        
        sc_ext_dm = ax_fert.scatter(x=xs, y=ys, s=s, c=cs, marker=marker, zorder=5)
        for x, y, e, c in zip(xs, ys, es, cs):
            eb = ax_fert.errorbar(x, y, yerr=e, lw=0,
                                  ecolor=c, elinewidth=1.0, capsize=1.5,
                                  marker=None, markersize=0.0,
                                  markeredgecolor='none',
                                  markeredgewidth=0.5,
                                  markerfacecolor=c,
                                  markerfacecoloralt=None, zorder=5)

        #ax_fert.axhline(0.75, color='#d62728', lw=1, zorder=6)
        if not last_page:
            ax_fert.set_xlabel('Core genes rank')
        ax_fert.set_ylabel('Fertility rate')
        ax_fert.set_yticks(np.linspace(0, 1, 5))
        ax_fert.set_xticks(xticks)
        ax_fert.set_xticklabels(xticklabels)
        ax_fert.set_ylim(-0.04, 1.04)
        ax_fert.set_xlim(xlimmin, xlimmax)
        ax_fert.grid(linewidth=0.5)
        ax_fert.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))


        if last_page:
            eb = ax_ctr.errorbar(range(0, len(dfc)), dfc['mean fert-rate'], yerr=dfc['std fert-rate'], lw=0,
                                 ecolor='#1f77b4', elinewidth=1.0, capsize=1.5,
                                 marker=marker, markersize=3.5,
                                 markeredgecolor='#1f77b4', markeredgewidth=0.5,
                                 markerfacecolor='#1f77b4', markerfacecoloralt=None, zorder=5)
            #ax_ctr.axhline(0.75, color='#d62728', lw=1, zorder=6)
            #ax_ctr.set_xlabel('Fertility Rate (Mean +/- SD)      ', fontsize='small', ha='center')
            ax_ctr.set_yticks(np.linspace(0, 1, 5))
            ax_ctr.set_yticklabels([])
            ax_ctr.set_xticks(range(0, len(dfc)))
            ax_ctr.set_xticklabels(['+ve C', '-ve C'], rotation=90, va='top', ha='center', fontsize='small')
            ax_ctr.set_ylim(-0.04, 1.04)
            ax_ctr.set_xlim(-.6, 1.6)
            ax_ctr.grid(linewidth=0.5)

            #
            # Legend
            #
            # Legend: FPKM
            norm_fpkm = mpl.colors.Normalize(vmin=0, vmax=11.79)
            cmap_fpkm = mpl.cm.Reds
            ax_fpkm_cb = fig.add_axes([0.67, 0.915, 0.16, 0.024])  # x, y, width, height
            ax_fpkm_cb.set_ylabel('FPKM', fontsize='small', ha='right') #, loc='left')
            cb_fpkm = mpl.colorbar.ColorbarBase(ax=ax_fpkm_cb, cmap=cmap_fpkm, norm=norm_fpkm, orientation='horizontal')
            cb_fpkm.set_label(r'log$_2$(FPKM+1)', fontsize='small')
            cb_fpkm.ax.tick_params(labelsize='x-small')

            # Legend: Phonetype (t = top; b = bottom)
            ax_code_cb_t = fig.add_axes([0.67, 0.50, 0.02, (0.038 * 5)])  # x, y, width, height
            ax_code_cb_b = fig.add_axes([0.83, 0.55, 0.02, (0.038 * 3)])  # x, y, width, height
            ax_code_cb_t.set_title('Phenotype code', fontsize='small', loc='left')
            cmap_code_t = mpl.colors.ListedColormap(list(code_color.values())[0:5][::-1])
            cmap_code_b = mpl.colors.ListedColormap(list(code_color.values())[5:8][::-1])
            bounds_code_t = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
            bounds_code_b = [5.0, 6.0, 7.0, 8.0]
            ticks_code_t = [0.5, 1.5, 2.5, 3.5, 4.5]
            ticks_code_b = [5.5, 6.5, 7.5]
            norm_code_t = mpl.colors.BoundaryNorm(bounds_code_t, cmap_code_t.N)
            norm_code_b = mpl.colors.BoundaryNorm(bounds_code_b, cmap_code_b.N)
            cb_code_t = mpl.colorbar.ColorbarBase(ax=ax_code_cb_t, cmap=cmap_code_t, norm=norm_code_t, boundaries=bounds_code_t, extend='neither', ticks=ticks_code_t, spacing='uniform', orientation='vertical')
            cb_code_b = mpl.colorbar.ColorbarBase(ax=ax_code_cb_b, cmap=cmap_code_b, norm=norm_code_b, boundaries=bounds_code_b, extend='neither', ticks=ticks_code_b, spacing='uniform', orientation='vertical')
            cb_code_t.set_ticklabels(list(code_label.values())[0:5][::-1])
            cb_code_b.set_ticklabels(list(code_label.values())[5:8][::-1])
            cb_code_t.ax.tick_params(labelsize='x-small')
            cb_code_b.ax.tick_params(labelsize='x-small')

            # Legend: Module
            ax_mod_cb = fig.add_axes([0.67, 0.05, 0.02, (0.038 * 9)])  # x, y, width, height
            ax_mod_cb.set_title('Module', fontsize='small', loc='left')
            cmap_mod = mpl.colors.ListedColormap(list(module_color.values())[::-1])
            bounds_mod = np.arange(0, len(module_color) + 1)
            ticks_mod = bounds_mod + 0.5
            norm_mod = mpl.colors.BoundaryNorm(bounds_mod, cmap_mod.N)
            cb_mod = mpl.colorbar.ColorbarBase(ax=ax_mod_cb, cmap=cmap_mod, norm=norm_mod, boundaries=bounds_mod, extend='neither', ticks=ticks_mod, spacing='uniform', orientation='vertical')
            cb_mod.set_ticklabels(list(module_label.values())[::-1])
            cb_mod.ax.tick_params(labelsize='x-small')


        plt.subplots_adjust(left=0.14, right=0.98, bottom=0.095, top=0.999, wspace=0.0, hspace=.9)
        file = 'images/img-core_DM_screened-{:d}.pdf'.format(page + 1)
        fig.savefig(file)

        #break
