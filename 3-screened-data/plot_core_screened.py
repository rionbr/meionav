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
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def value_to_color(x, cmap, norm):
    rgb = cmap(norm(x))[:3]
    return colors.rgb2hex(rgb)


if __name__ == '__main__':

    # Load genes
    df = pd.read_csv('../2-core_genes/results/pipeline-core/DM_meiotic_genes.csv', index_col=0, usecols=['id_gene', 'gene'])

    # Load Screened data
    dfs = pd.read_csv('data/core_DM_screened_2020-01-07.csv', index_col=0)

    # Load FPKM data
    dfFPKM1 = pd.read_csv('../1-diff-gene-exp/data/DM/DM_Spermatocytes_FPKM_sample1.csv', index_col=0, usecols=['Gene ID', 'FPKM'])
    dfFPKM2 = pd.read_csv('../1-diff-gene-exp/data/DM/DM_Spermatocytes_FPKM_sample2.csv', index_col=0, usecols=['Gene ID', 'FPKM'])

    dfs_only = dfs.loc[~dfs['FT1 eggs'].isnull(), :]

    status_cats = ['Screened', 'To be crossed', 'Pending', 'Reorder']
    dfs['Status'] = pd.Categorical(dfs['Status'], categories=status_cats, ordered=True)
    df['Status'] = dfs['Status']

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

    print(dfs.head())
    print(dfs.loc[dfs['MM pheno code'].str.len() > 1, :])
    df['RNAi'] = dfs['Previous ref to RNAi working?'].apply(lambda x: 'Yes' if 'Yes' in x else 'No')
    df['our-DM-code'] = dfs['Our DM pheno code']
    df['ext-DM-code'] = dfs['Others DM pheno code']
    df['ext-MM-code'] = dfs['MM pheno code']
    df['ext-HS-code'] = dfs['HS pheno code']
    print(df.head)

    # FPKM
    df['FPKM1'] = dfFPKM1['FPKM']
    df['FPKM2'] = dfFPKM2['FPKM']
    df['FPKM'] = df[['FPKM1', 'FPKM2']].mean(axis='columns')
    df['logFPKM'] = df['FPKM'].apply(lambda x: np.log2(x + 1))

    df = df.sort_values(['Status', 'mean fert-rate', 'gene'], ascending=[True, True, True]).reset_index()

    maxfpkm, minfpkm = df['logFPKM'].max(), df['logFPKM'].min()

    codemap = {

    }
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
        'A': '#1f77b4',
        'B': '#ff7f0e',
        'C': '#2ca02c',
        'D': '#d62728',
        'E': '#9467bd',
        'F': '#8c564b',
        'G': '#e377c2',
        'H': '#7f7f7f'
    }

    # print(df.head())
    # print(df.tail())
    print("logFPKM: {:.2f}/{:.2f}".format(minfpkm, maxfpkm))

    #
    # Plot Page 1
    #
    print("Plotting")
    n_per_page = 150
    dfg = df.groupby(np.arange(len(df)) // n_per_page)
    number_of_pages = len(dfg)
    for page, dft in dfg:
        dft = dft.sort_values(['Status', 'mean fert-rate', 'gene'], ascending=[False, False, False]).reset_index()
        n_this_page = len(dft)
        print("> Page: {:d} of {:d}".format(page + 1, number_of_pages))
        print("> Points in this page: {:d}".format(n_this_page))

        fig = plt.figure(figsize=(4, 11))
        # fig.suptitle('Core metazoan meiotic genes'.format(page, number_of_pages))

        gs = gridspec.GridSpec(nrows=n_per_page, ncols=14)
        ax_fert = plt.subplot(gs[:n_this_page, 0:8])
        ax_fpkm = plt.subplot(gs[:n_this_page, 8])
        ax_rnai = plt.subplot(gs[:n_this_page, 9])
        ax_our_dm = plt.subplot(gs[:n_this_page, 10])
        ax_ext_dm = plt.subplot(gs[:n_this_page, 11])
        ax_ext_mm = plt.subplot(gs[:n_this_page, 12])
        ax_ext_hs = plt.subplot(gs[:n_this_page, 13])

        adjustable = 'datalim'
        aspect = 'auto'
        ax_fert.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_fpkm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_rnai.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_our_dm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_ext_dm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_ext_mm.set(adjustable=adjustable, aspect=aspect, anchor='NE')
        ax_ext_hs.set(adjustable=adjustable, aspect=aspect, anchor='NE')

        norm_fpkm = mpl.colors.Normalize(vmin=minfpkm, vmax=maxfpkm)
        cmap_fpkm = mpl.cm.Reds

        rotation = 50
        s = 12
        marker = 's'
        #
        # DM Expression Values
        #
        eb = ax_fert.errorbar(dft['mean fert-rate'], range(0, len(dft)), xerr=dft['std fert-rate'], lw=0,
                              ecolor='#3182bd', elinewidth=1.0, capsize=2,
                              marker=marker, markersize=3.5,
                              markeredgecolor='#3182bd', markeredgewidth=0.5,
                              markerfacecolor='#6baed6', markerfacecoloralt=None, zorder=5)
        ax_fert.axvline(0.75, color='#d62728', lw=1, zorder=6)
        ax_fert.set_xlabel('Fertility Rate (Mean +/- SD)      ', fontsize='small', ha='center')
        ax_fert.set_xticks(np.linspace(0, 1, 5))
        ax_fert.set_xticklabels(np.linspace(0, 1, 5), fontsize='small', rotation=0)
        ax_fert.set_yticks(range(0, len(dft)))
        ax_fert.set_yticklabels(dft['gene'], rotation=0, va='center', ha='right', fontsize='xx-small')
        ax_fert.set_xlim(-0.04, 1.04)
        ax_fert.set_ylim(-1, len(dft))
        ax_fert.grid(linewidth=0.5)

        #
        # Expression (FPKM)
        #
        y = dft['logFPKM'].index
        x = np.zeros(len(y))
        c = dft['logFPKM'].apply(value_to_color, args=(cmap_fpkm, norm_fpkm))

        sc_fpkm = ax_fpkm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_fpkm.set_xticks([0])
        ax_fpkm.set_xticklabels(['FPKM'], fontsize='small', rotation=rotation, rotation_mode='anchor', va='top', ha='right')
        ax_fpkm.set_yticks(range(0, len(dft)))
        ax_fpkm.tick_params(axis='y', which='major', length=1.5)
        ax_fpkm.set_yticklabels([])
        ax_fpkm.set_xlim(-0.2, 0.2)  # Adjusting this makes the plot shrink
        ax_fpkm.set_ylim(-1, len(dft))

        #
        # RNAi
        #
        data_rnai = dft.loc[dft['RNAi'] == 'Yes', 'RNAi']
        y = data_rnai.index
        x = np.zeros(len(y))
        c = '#17becf'

        sc_rnai = ax_rnai.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_rnai.set_xticks([0])
        ax_rnai.set_xticklabels(['RNAi efficiency'], fontsize='small', rotation=rotation, rotation_mode='anchor', va='top', ha='right')
        ax_rnai.set_yticks(range(0, len(dft)))
        ax_rnai.tick_params(axis='y', which='major', length=1.5)
        ax_rnai.set_yticklabels([])
        ax_rnai.set_xlim(-0.2, 0.2)
        ax_rnai.set_ylim(-1, len(dft))
        # ax_rnai.grid(axis='y', linewidth=0.5)

        #
        # Our DM Phenotype
        #
        data_our_dm = dft.loc[~dft['our-DM-code'].isnull(), 'our-DM-code']
        y = data_our_dm.index
        x = np.zeros(len(y))
        c = data_our_dm.map(code_color)

        sc_our_dm = ax_our_dm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_our_dm.set_xticks([0])
        ax_our_dm.set_xticklabels(['*DM phenotype'], fontsize='small', rotation=rotation, rotation_mode='anchor', va='top', ha='right')
        ax_our_dm.set_yticks(range(0, len(dft)))
        ax_our_dm.tick_params(axis='y', which='major', length=1.5)
        ax_our_dm.set_yticklabels([])
        ax_our_dm.set_xlim(-0.5, 0.5)
        ax_our_dm.set_ylim(-1, len(dft))
        # ax_our_dm.grid(axis='y', linewidth=0.5)

        #
        # External DM Phenotype
        #
        data_ext_dm = dft.loc[~dft['ext-DM-code'].isnull(), 'ext-DM-code']
        y = data_ext_dm.index
        x = np.zeros(len(y))
        c = data_ext_dm.map(code_color)

        sc_ext_dm = ax_ext_dm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_ext_dm.set_xticks([0])
        ax_ext_dm.set_xticklabels(['DM phenotype'], fontsize='small', rotation=rotation, rotation_mode='anchor', va='top', ha='right')
        ax_ext_dm.set_yticks(range(0, len(dft)))
        ax_ext_dm.tick_params(axis='y', which='major', length=1.5)
        ax_ext_dm.set_yticklabels([])
        ax_ext_dm.set_xlim(-0.5, 0.5)
        ax_ext_dm.set_ylim(-1, len(dft))
        # ax_ext_dm.grid(axis='y', linewidth=0.5)

        #
        # External MM Phenotype
        #
        # (these lines solve the problem when an MM phenotype has two codes, e.g., A/B)
        data_ext_mm = dft.loc[~dft['ext-MM-code'].isnull(), 'ext-MM-code']
        if len(data_ext_mm):
            # print(data_ext_mm)
            data_tmp = data_ext_mm.str.split('/').apply(pd.Series)
            # print(data_tmp)
            data_tmp = pd.melt(data_tmp.reset_index(), id_vars='index', value_vars=data_tmp.columns.tolist())
            # print(data_tmp)
            data_tmp = data_tmp.set_index('index').dropna(subset=['value'])
            # print(data_tmp)
            data_tmp.loc[(data_tmp.index.duplicated(keep=False) & (data_tmp['variable'] == 0)), 'variable'] = -0.2
            data_tmp.loc[(data_tmp.index.duplicated(keep=False) & (data_tmp['variable'] == 1)), 'variable'] = +0.2
            #
            y = data_tmp.index.values
            x = data_tmp['variable'].values
            c = data_tmp['value'].map(code_color).values
            sc_ext_mm = ax_ext_mm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_ext_mm.set_xticks([0])
        ax_ext_mm.set_xticklabels(['MM phenotype'], fontsize='small', rotation=rotation, rotation_mode='anchor', va='top', ha='right')
        ax_ext_mm.set_yticks(range(0, len(dft)))
        ax_ext_mm.tick_params(axis='y', which='major', length=1.5)
        ax_ext_mm.set_yticklabels([])
        ax_ext_mm.set_xlim(-0.5, 0.5)
        ax_ext_mm.set_ylim(-1, len(dft))
        # ax_ext_mm.grid(axis='y', linewidth=0.5)

        #
        # External HS Phenotype
        #
        data_ext_hs = dft.loc[~dft['ext-HS-code'].isnull(), 'ext-HS-code']
        y = data_ext_hs.index
        x = np.zeros(len(y))
        c = data_ext_hs.map(code_color)

        sc_ext_hs = ax_ext_dm.scatter(x=x, y=y, s=s, c=c, marker=marker, zorder=5)
        ax_ext_hs.set_xticks([0])
        ax_ext_hs.set_xticklabels(['HS phenotype'], fontsize='small', rotation=rotation, rotation_mode='anchor', va='top', ha='right')
        ax_ext_hs.set_yticks(range(0, len(dft)))
        ax_ext_hs.tick_params(axis='y', which='major', length=1.5)
        ax_ext_hs.set_yticklabels([])
        ax_ext_hs.set_xlim(-0.5, 0.5)
        ax_ext_hs.set_ylim(-1, len(dft))
        # ax_ext_hs.grid(axis='y', linewidth=0.5)

        plt.subplots_adjust(left=0.2, right=0.96, bottom=0.08, top=0.96, wspace=0.2, hspace=0.0)
        file = 'images/img-core_DM_screened-{:d}.pdf'.format(page + 1)
        fig.savefig(file)
