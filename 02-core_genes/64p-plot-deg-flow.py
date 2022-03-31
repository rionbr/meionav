# coding=utf-8
# Author: Rion B Correia
# Date: Nov 05, 2020
#
# Description: Plots an overview of the FPKM across the three species
#
#
import math
import numpy as np
import pandas as pd
idx = pd.IndexSlice
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
#
mpl.rc('font', size=10)  # controls default text sizes
mpl.rc('axes', titlesize=12)  # fontsize of the axes title
mpl.rc('axes', labelsize=12)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=10)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=10)  # fontsize of the tick labels
mpl.rc('legend', fontsize=10)  # legend fontsize
mpl.rc('figure', titlesize=12)  # fontsize of the figure title
#
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
from utils import ensurePathExists
from sklearn import preprocessing
import ast
from matplotlib.patches import ConnectionPatch


def load_species_de_df(species='DM'):
    rDEFile = 'results/DE/{species:s}-DE_genes.csv.gz'.format(species=species)
    #
    df = pd.read_csv(rDEFile, index_col=0)
    #
    # Entry
    #
    if species == 'DM':
        col_entry = 'Middle_vs_Apical'
        col_var_entry = 'MiddleApical'
    else:
        col_entry = 'Cyte_vs_Gonia'
        col_var_entry = 'CyteGonia'

    col_logFC_entry = 'logFC_{col:s}'.format(col=col_var_entry)
    col_FDR_entry = 'FDR_{col:s}'.format(col=col_var_entry)
    #
    df_entry = df.loc[
        (df[col_entry] == True) &
        ((df[col_logFC_entry] >= 1) | (df[col_logFC_entry] <= -1)) &
        (df[col_FDR_entry] <= 0.05)
        , :].copy()
    #
    df_entry['rank'] = df_entry[col_logFC_entry].rank(method='dense', ascending=False)
    df_entry.sort_values('rank', ascending=True, inplace=True)

    #
    # Exit
    #
    if species == 'DM':
        col_exit = 'Basal_vs_Middle'
        col_var_exit = 'BasalMiddle'
    else:
        col_exit = 'Tid_vs_Cyte'
        col_var_exit = 'TidCyte'

    col_logFC_exit = 'logFC_{col:s}'.format(col=col_var_exit)
    col_FDR_exit = 'FDR_{col:s}'.format(col=col_var_exit)

    # Only DGE (up or down) genes
    df_exit = df.loc[
        (df[col_exit] == True) &
        ((df[col_logFC_exit] >= 1) | (df[col_logFC_exit] <= -1)) &
        (df[col_FDR_exit] <= 0.05)
        , :].copy()
    #
    df_exit['rank'] = df_exit[col_logFC_exit].rank(method='dense', ascending=False)
    df_exit.sort_values('rank', ascending=True, inplace=True)
    #
    return df_entry, df_exit


def load_species_core_df(species='DM'):
    rCoreFile = 'results/pipeline-core/{species:s}_meiotic_genes.csv'.format(species=species)
    dfc = pd.read_csv(rCoreFile, index_col=0)

    if species == 'DM':
        col_entry = 'Middle_vs_Apical'
        col_var_entry = 'MiddleApical'
        col_exit = 'Basal_vs_Middle'
        col_var_exit = 'BasalMiddle'
    else:
        col_entry = 'Cyte_vs_Gonia'
        col_var_entry = 'CyteGonia'
        col_exit = 'Tid_vs_Cyte'
        col_var_exit = 'TidCyte'

    col_logFC_entry = 'logFC_{col:s}'.format(col=col_var_entry)
    col_FDR_entry = 'FDR_{col:s}'.format(col=col_var_entry)
    #
    col_logFC_exit = 'logFC_{col:s}'.format(col=col_var_exit)
    col_FDR_exit = 'FDR_{col:s}'.format(col=col_var_exit)

    if species == 'DM':
        dfc_entry = dfc.loc[(dfc[col_entry] == True), :].copy()
        dfc_exit = dfc.loc[(dfc[col_exit] == True), :].copy()
    else:
        dfc_entry = dfc.loc[(dfc[col_entry] == True) & (dfc[col_logFC_entry] >= 1) & (dfc[col_FDR_entry] <= 0.05), :].copy()
        
        dfc_exit = dfc.loc[(dfc[col_exit] == True) & (dfc[col_logFC_exit] <= -1) & (dfc[col_FDR_exit] <= 0.05), :].copy()

    return dfc_entry, dfc_exit


def load_species_fpkm_dfs(species, df_entry, df_exit):
    ldfs = []
    #
    for celltype in ['spermatogonia', 'spermatocyte', 'spermatid']:

        rCSVFile = 'results/FPKM/{species:s}/{species:s}-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype, species=species)
        df = pd.read_csv(rCSVFile, index_col=0)

        df['id_string'] = df['id_string'].apply(lambda x: ast.literal_eval(x) if str(x).startswith('[') else [x])
        if celltype == 'spermatocyte':
            dict_id_gene_id_string = df.explode('id_string').dropna().reset_index().set_index('id_string')['id_gene'].to_dict()

        # Drop Duplicated
        df = df.loc[~df.index.duplicated(), :]
        # Calc Log(FPKM)
        df['logFPKM'] = df['FPKM'].apply(lambda x: np.log2(x + 1))
        # Index to MultiIndex
        index = pd.MultiIndex.from_arrays([[celltype] * 6, df.columns])
        df.columns = index

        ldfs.append(df)
    dfc = pd.concat(ldfs, axis='columns', verify_integrity=True).fillna(0)

    dfc = dfc.loc[:, idx[:, 'logFPKM']]
    dfc.columns = dfc.columns.droplevel(level=1)

    #
    # Entry
    #
    # index
    dfc_entry = dfc.loc[(dfc.index.isin(df_entry.index)), ['spermatogonia', 'spermatocyte']]
    dfcn_entry = dfc_entry
    dfcn_entry = pd.DataFrame(preprocessing.normalize(dfc_entry.values, norm='l2'), columns=dfc_entry.columns, index=dfc_entry.index)
    # rank/sort/y
    dfcn_entry['rank'] = dfcn_entry.index.map(df_entry['rank'].to_dict())
    dfcn_entry.sort_values('rank', ascending=True, inplace=True)
    dfcn_entry['y'] = np.arange(len(dfcn_entry))

    #
    # Exit
    #
    dfc_exit = dfc.loc[(dfc.index.isin(df_exit.index)), ['spermatocyte', 'spermatid']]
    dfcn_exit = dfc_exit
    dfcn_exit = pd.DataFrame(preprocessing.normalize(dfc_exit.values, norm='l2'), columns=dfc_exit.columns, index=dfc_exit.index)
    # rank/sort/y
    dfcn_exit['rank'] = dfcn_exit.index.map(df_exit['rank'].to_dict())
    dfcn_exit.sort_values('rank', ascending=True, inplace=True)
    dfcn_exit['y'] = np.arange(len(dfcn_exit))

    #
    return dfcn_entry, dfcn_exit, dict_id_gene_id_string


def load_meta_genes():
    df = pd.read_csv('results/meta-genes/meta-spermatocyte-genes.csv.gz')
    df['id_string_HS'] = df['id_string_HS'].str.split(',')
    df['id_string_MM'] = df['id_string_MM'].str.split(',')
    df['id_string_DM'] = df['id_string_DM'].str.split(',')

    # Explode
    df = df.explode('id_string_HS').explode('id_string_MM').explode('id_string_DM')

    df['id_gene_HS'] = df['id_string_HS'].map(dict_hs_gene_string)
    df['id_gene_MM'] = df['id_string_MM'].map(dict_mm_gene_string)
    df['id_gene_DM'] = df['id_string_DM'].map(dict_dm_gene_string)

    return df


if __name__ == '__main__':

    #
    # HS
    #

    # DE
    df_de_hs_entry, df_de_hs_exit = load_species_de_df(species='HS')
    # FPKM
    df_fpkm_hs_entry, df_fpkm_hs_exit, dict_hs_gene_string = load_species_fpkm_dfs(species='HS', df_entry=df_de_hs_entry, df_exit=df_de_hs_exit)

    # Invert downregulated fpkm
    df_fpkm_hs_entry.loc[(df_de_hs_entry['logFC_CyteGonia'] < -1), ['spermatogonia', 'spermatocyte']] *= -1
    df_fpkm_hs_exit.loc[(df_de_hs_exit['logFC_TidCyte'] < -1), ['spermatocyte', 'spermatid']] *= -1

    # Core
    df_core_hs_entry, df_core_hs_exit = load_species_core_df(species='HS')

    #
    # MM
    #

    # DE
    df_de_mm_entry, df_de_mm_exit = load_species_de_df(species='MM')
    # FPKM
    df_fpkm_mm_entry, df_fpkm_mm_exit, dict_mm_gene_string = load_species_fpkm_dfs(species='MM', df_entry=df_de_mm_entry, df_exit=df_de_mm_exit)

    # Invert downregulated fpkm
    df_fpkm_mm_entry.loc[(df_de_mm_entry['logFC_CyteGonia'] < -1), ['spermatogonia', 'spermatocyte']] *= -1
    df_fpkm_mm_exit.loc[(df_de_mm_exit['logFC_TidCyte'] < -1), ['spermatocyte', 'spermatid']] *= -1

    # Core
    df_core_mm_entry, df_core_mm_exit = load_species_core_df(species='MM')

    #
    # DM
    #

    # DE
    df_de_dm_entry, df_de_dm_exit = load_species_de_df(species='DM')
    # FPKM
    df_fpkm_dm_entry, df_fpkm_dm_exit, dict_dm_gene_string = load_species_fpkm_dfs(species='DM', df_entry=df_de_dm_entry, df_exit=df_de_dm_exit)

    # Invert downregulated fpkm
    df_fpkm_dm_entry.loc[(df_de_dm_entry['logFC_MiddleApical'] < -1), ['spermatogonia', 'spermatocyte']] *= -1
    df_fpkm_dm_exit.loc[(df_de_dm_exit['logFC_BasalMiddle'] < -1), ['spermatocyte', 'spermatid']] *= -1

    # Core
    df_core_dm_entry, df_core_dm_exit = load_species_core_df(species='DM')

    # Load Meta Genes
    df_meta = load_meta_genes()

    # Subselect to Plot
    col_hs = ['id_gene_HS', 'id_string_HS']
    col_mm = ['id_gene_MM', 'id_string_MM']
    col_dm = ['id_gene_DM', 'id_string_DM']

    # HS-MM
    df_hs_mm_entry = df_meta.loc[(df_meta['id_gene_HS'].isin(df_core_hs_entry.index)) & (df_meta['id_gene_MM'].isin(df_core_mm_entry.index)), col_hs + col_mm].drop_duplicates()
    df_hs_mm_exit = df_meta.loc[(df_meta['id_gene_HS'].isin(df_core_hs_exit.index)) & (df_meta['id_gene_MM'].isin(df_core_mm_exit.index)), col_hs + col_mm].drop_duplicates()
    #
    df_hs_mm_entry['y-HS'] = df_hs_mm_entry['id_gene_HS'].map(df_fpkm_hs_entry['y'])
    df_hs_mm_entry['y-MM'] = df_hs_mm_entry['id_gene_MM'].map(df_fpkm_mm_entry['y'])
    df_hs_mm_exit['y-HS'] = df_hs_mm_exit['id_gene_HS'].map(df_fpkm_hs_exit['y'])
    df_hs_mm_exit['y-MM'] = df_hs_mm_exit['id_gene_MM'].map(df_fpkm_mm_exit['y'])

    # MM-DM
    df_mm_dm_entry = df_meta.loc[df_meta['id_gene_MM'].isin(df_core_mm_entry.index) & df_meta['id_gene_DM'].isin(df_de_dm_entry.index.tolist() + df_de_dm_exit.index.to_list()), col_mm + col_dm].drop_duplicates()
    df_mm_dm_exit = df_meta.loc[df_meta['id_gene_MM'].isin(df_core_mm_exit.index) & df_meta['id_gene_DM'].isin(df_de_dm_entry.index.tolist() + df_de_dm_exit.index.to_list()), col_mm + col_dm].drop_duplicates()
    #
    df_mm_dm_entry['y-MM'] = df_mm_dm_entry['id_gene_MM'].map(df_fpkm_mm_entry['y'])
    df_mm_dm_entry['y-DM'] = df_mm_dm_entry['id_gene_DM'].map(df_fpkm_dm_entry['y'])
    df_mm_dm_exit['y-MM'] = df_mm_dm_exit['id_gene_MM'].map(df_fpkm_mm_exit['y'])
    df_mm_dm_exit['y-DM'] = df_mm_dm_exit['id_gene_DM'].map(df_fpkm_dm_exit['y'])

    # DEBUG
    # display(df_hs_mm_entry.loc[df_hs_mm_entry['id_gene_HS'].isin(df_core_hs_exit.index)])
    # print('df_hs_mm_entry:',df_hs_mm_entry.shape)
    # print('df_hs_mm_exit:',df_hs_mm_exit.shape)
    # print('df_mm_dm_entry:',df_mm_dm_entry.shape)
    # print('df_mm_dm_exit:',df_mm_dm_exit.shape)
    aspect = 1/900
    #
    # Plot Entry
    #
    fig = plt.figure(figsize=(3.7, 2.9))
    #
    gs = GridSpec(nrows=100, ncols=3)
    ax_hs = plt.subplot(gs[0:46, 0], anchor='N')
    ax_mm = plt.subplot(gs[0:, 1], anchor='N')
    ax_dm = plt.subplot(gs[0:61, 2], anchor='N')
    #
    axes = [ax_hs, ax_mm, ax_dm]
    #
    #caxu = plt.axes([0.09, (0.12), 0.18, 0.018])
    #caxd = plt.axes([0.09, (0.12 - 0.020), 0.18, 0.018])
    #caxl = plt.axes([0.73, (0.12), 0.18, 0.018])

    species_label = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    # red = '#d62728'
    # blue = '#1f77b4'
    # green = '#2ca02c'
    cmap = 'coolwarm'
    vmin = -1
    vmax = 1
    norm = Normalize(vmin=vmin, vmax=vmax)

    dfs = [df_fpkm_hs_entry, df_fpkm_mm_entry, df_fpkm_dm_entry]
    species = ['HS', 'MM', 'DM']
    
    for ax, df, species in zip(axes, dfs, species):
        values = df[['spermatogonia', 'spermatocyte']].values
        im = ax.imshow(values, cmap=cmap, norm=norm, aspect=aspect, interpolation='nearest')

        yticks = np.arange(0, len(df), 1000)

        ax.set_title(species_label[species])
        ax.set_xticks([0, 1])
        ax.set_yticks(yticks)
        ax.get_yaxis().set_major_formatter(lambda x, p: format(int(x), ','))
        ax.set_xticklabels(['Sg.', 'Sc.'])

    for i, row in df_hs_mm_entry.iterrows():
        y_hs = row['y-HS']
        y_mm = row['y-MM']
        con = ConnectionPatch((1.5, y_hs), (-0.5, y_mm), coordsA='data', coordsB='data', axesA=ax_hs, axesB=ax_mm, color='#2ca02c', alpha=.3, zorder=-1)
        fig.add_artist(con)

    for i, row in df_mm_dm_entry.iterrows():
        y_mm = row['y-MM']
        y_dm = row['y-DM']
        con = ConnectionPatch((1.5, y_mm), (-0.5, y_dm), coordsA='data', coordsB='data', axesA=ax_mm, axesB=ax_dm, color='#2ca02c', alpha=.3, zorder=-1)
        fig.add_artist(con)
    
    """
    cbu = plt.colorbar(mpl.cm.ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='Reds'), ticks=[], orientation='horizontal', cax=caxu)
    cbd = plt.colorbar(mpl.cm.ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='Blues'), ticks=[0.0, 0.5, 1.0], orientation='horizontal', cax=caxd)
    #caxu.yaxis.tick_left()
    caxu.set_title('norm[log(FPKM+1)]', fontsize=10)
    #
    caxl.plot((0, 1), (0, 0), color='#2ca02c')
    caxl.set_title('Orthologs', fontsize=10)
    #caxl.xaxis.set_label_position("right")
    #caxl.yaxis.tick_right()
    caxl.tick_params(axis='both', which='both', bottom=False, left=False, right=False, labelleft=False, labelbottom=False, labelright=False)
    """
    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.10, top=0.92, wspace=1.2)

    wIMGfile = 'images/deg-flow/img-deg-flow-entry.pdf'
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=300)
    plt.close()

    #
    # Plot Exit
    #
    fig = plt.figure(figsize=(3.7, 2.9))
    #
    gs = GridSpec(nrows=100, ncols=3)
    ax_hs = plt.subplot(gs[0:88, 0], anchor='N')
    ax_mm = plt.subplot(gs[0:, 1], anchor='N')
    ax_dm = plt.subplot(gs[0:19, 2], anchor='N')
    #
    axes = [ax_hs, ax_mm, ax_dm]
    #
    caxl = plt.axes([0.75, (0.27), 0.18, 0.018])
    caxu = plt.axes([0.75, (0.27 - 0.13), 0.18, 0.018])
    caxd = plt.axes([0.75, (0.27 - 0.13 - 0.020), 0.18, 0.018])

    species_label = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    # red = '#d62728'
    # blue = '#1f77b4'
    # green = '#2ca02c'
    cmap = 'coolwarm'  # LinearSegmentedColormap.from_list(name='custom-cmap', colors=['white', '#d62728'])
    vmin = -1
    vmax = 1
    norm = Normalize(vmin=vmin, vmax=vmax)

    dfs = [df_fpkm_hs_exit, df_fpkm_mm_exit, df_fpkm_dm_exit]
    species = ['HS', 'MM', 'DM']
    for ax, df, species in zip(axes, dfs, species):
        values = df[['spermatocyte', 'spermatid']].values
        im = ax.imshow(values, cmap=cmap, norm=norm, aspect=aspect, interpolation='nearest')

        yticks = np.arange(0, len(df), 1000)

        ax.set_title(species_label[species])
        ax.set_xticks([0, 1])
        ax.set_yticks(yticks)
        ax.get_yaxis().set_major_formatter(lambda x, p: format(int(x), ','))
        ax.set_xticklabels(['Sc.', 'St.'])

    for i, row in df_hs_mm_exit.iterrows():
        y_hs = row['y-HS']
        y_mm = row['y-MM']
        con = ConnectionPatch((1.5, y_hs), (-0.5, y_mm), coordsA='data', coordsB='data', axesA=ax_hs, axesB=ax_mm, color='#2ca02c', alpha=.3, zorder=-1)
        fig.add_artist(con)

    for i, row in df_mm_dm_exit.iterrows():
        y_mm = row['y-MM']
        y_dm = row['y-DM']
        con = ConnectionPatch((1.5, y_mm), (-0.5, y_dm), coordsA='data', coordsB='data', axesA=ax_mm, axesB=ax_dm, color='#2ca02c', alpha=.3, zorder=-1)
        fig.add_artist(con)

    #ax_dm.yaxis.tick_right()
    cbu = plt.colorbar(mpl.cm.ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='Reds'), ticks=[], orientation='horizontal', cax=caxu)
    cbd = plt.colorbar(mpl.cm.ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='Blues'), ticks=[0.0, 0.5, 1.0], orientation='horizontal', cax=caxd)
    #caxu.yaxis.tick_left()
    caxu.set_title('Exp. level', fontsize=10)
    #
    caxl.axhline(y=0.5, color='#2ca02c', lw=1.5)
    caxl.axis('off')
    caxl.set_title('Orthologs', fontsize=10, pad=4)
    caxl.tick_params(axis='both', which='both', bottom=False, left=False, right=False, labelleft=False, labelbottom=False, labelright=False)

    plt.subplots_adjust(left=0.14, right=0.96, bottom=0.10, top=0.92, wspace=1.2)

    wIMGfile = 'images/deg-flow/img-deg-flow-exit.pdf'
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=300)
    plt.close()