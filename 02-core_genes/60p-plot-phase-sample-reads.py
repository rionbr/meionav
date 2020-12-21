# coding=utf-8
# Author: Rion B Correia
# Date: Nov 05, 2020
#
# Description: Plots an overview of the FPKM across the three species
#
#
import math
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
from utils import ensurePathExists
from sklearn import preprocessing


phase_DM_translate = {
    'Cyte_vs_Gonia': 'Middle_vs_Apical',  # entry
    'Tid_vs_Cyte': 'Basal_vs_Middle'  # exit
}


def get_diff_expr_genes(specie, phase, updown):
    if phase not in ['Cyte_vs_Gonia', 'Tid_vs_Cyte']:
        raise TypeError('Input not valid')
    #
    if (specie == 'DM'):
        phaset = phase_DM_translate.get(phase)
    else:
        phaset = phase
    #
    minLogFC = math.log2(2)
    maxFDR = 0.05
    #
    file = '../01-diff-gene-exp/results/{specie:s}/{specie:s}-DGE_{phase:s}.csv'.format(specie=specie, phase=phaset)
    df = pd.read_csv(file, index_col=0)
    #
    dfU = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] >= 0), :]
    dfD = df.loc[(df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC) & (df['logFC'] <= 0), :]
    #
    if updown == 'Up':
        return set(dfU.index.tolist())
    elif updown == 'Down':
        return set(dfD.index.tolist())
    elif updown == 'Both':
        return set(dfU.index.tolist() + dfD.index.tolist())
    else:
        raise Exception('Up or Down?')


def probnorm(df):
    dfp = df.copy()
    #
    dfp = dfp.apply(lambda x: x / x.sum(), axis='index')
    #
    dfp = pd.DataFrame(preprocessing.normalize(dfp, norm='l2'), columns=dfp.columns, index=dfp.index)
    return dfp


def plot(df, species, phase, figsize=(4, 3)):
    fig, ax = plt.subplots(figsize=figsize)
    #
    speciest = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}[species]
    color = {'HS': '#2ca02c', 'MM': '#7f7f7f', 'DM': '#ff7f0e'}[species]
    positions = list(range(len(df.columns)))
    #
    ax.plot(positions, df.sample(frac=0.3).T, color=color, lw=1, alpha=0.05, zorder=5, rasterized=False)
    ax.boxplot(df.values, positions=positions, meanline=True, notch=False, widths=0.6, showfliers=False, medianprops={'color': 'black', 'lw': 1}, zorder=4)
    #
    ax.set_title('{species:s} meiotic {phase:s}'.format(species=speciest, phase=phase))
    ax.set_xticks(positions)
    ax.set_xticklabels(df.columns, rotation=90)
    ax.set_ylabel('Read counts (L2 norm)')
    #ax.set_xlabel('Sample')

    plt.tight_layout()
    #plt.subplots_adjust(left=0.17, bottom=0.17, right=0.97, top=0.90)
    plt.subplots_adjust(bottom=0.32, top=0.90)
    wIMGfile = 'images/phase-reads/img-{species:s}-meiotic-{phase:s}.pdf'.format(species=species, phase=phase)
    ensurePathExists(wIMGfile)
    fig.savefig(wIMGfile, dpi=150)


if __name__ == '__main__':

    print('Plotting Human')
    hs_entry_deg = get_diff_expr_genes('HS', 'Cyte_vs_Gonia', 'Up')
    hs_exit_deg = get_diff_expr_genes('HS', 'Tid_vs_Cyte', 'Down')

    df_hs = pd.read_csv('../01-diff-gene-exp/data/DGE/HS/HS_RawCounts_AllSperm.csv', index_col=0)
    #
    dict_cols_gonia = {
        'Spermatogonia_Ad_1': 'Gonia Ad-1',
        'Spermatogonia_Ad_2': 'Gonia Ad-2',
        'Spermatogonia_Ad_3': 'Gonia Ad-3',
        'Spermatogonia_Ad_4': 'Gonia Ad-4',
        'Spermatogonia_Ad_5': 'Gonia Ad-5',
        'Spermatogonia_Ad_6': 'Gonia Ad-6',
        'Spermatogonia_Ap_1': 'Gonia Ap-1',
        'Spermatogonia_Ap_2': 'Gonia Ap-2',
        'Spermatogonia_Ap_3': 'Gonia Ap-3',
        'Spermatogonia_Ap_4': 'Gonia Ap-4',
        'Spermatogonia_Ap_5': 'Gonia Ap-5',
    }
    dict_cols_cyte = {
        'Spermatocytes_Early_1': 'Cyte E-1',
        'Spermatocytes_Early_2': 'Cyte E-2',
        'Spermatocytes_Early_3': 'Cyte E-3',
        'Spermatocytes_Early_4': 'Cyte E-4',
        'Spermatocytes_Early_5': 'Cyte E-5',
        'Spermatocytes_Early_6': 'Cyte E-6',
        'Spermatocytes_Late_1': 'Cyte L-1',
        'Spermatocytes_Late_2': 'Cyte L-2',
        'Spermatocytes_Late_3': 'Cyte L-3',
        'Spermatocytes_Late_4': 'Cyte L-4',
        'Spermatocytes_Late_5': 'Cyte L-5',
        'Spermatocytes_Late_6': 'Cyte L-6',
    }
    dict_cols_tid = {
        'Spermatids_1': 'Tid-1',
        'Spermatids_2': 'Tid-2',
        'Spermatids_3': 'Tid-3',
        'Spermatids_4': 'Tid-4',
        'Spermatids_5': 'Tid-5',
        'Spermatids_6': 'Tid-6',
    }
    dict_all_cols = {**dict_cols_gonia, **dict_cols_cyte, **dict_cols_tid}
    cols_gonia = list(dict_cols_gonia.values())
    cols_cyte = list(dict_cols_cyte.values())
    cols_tid = list(dict_cols_tid.values())
    df_hs = df_hs.rename(columns=dict_all_cols)
    #
    df_hs_entry = df_hs.loc[df_hs.index.isin(hs_entry_deg), cols_gonia + cols_cyte]
    df_hs_exit = df_hs.loc[df_hs.index.isin(hs_exit_deg), cols_cyte + cols_tid]
    #
    #df_hs_entry = pd.DataFrame(preprocessing.normalize(df_hs_entry, norm='l2'), columns=df_hs_entry.columns, index=df_hs_entry.index)
    #df_hs_exit = pd.DataFrame(preprocessing.normalize(df_hs_exit, norm='l2'), columns=df_hs_exit.columns, index=df_hs_exit.index)
    df_hs_entry = probnorm(df_hs_entry)
    df_hs_exit = probnorm(df_hs_exit)
    #
    plot(df_hs_entry, 'HS', 'entry', figsize=(6, 3))
    plot(df_hs_exit, 'HS', 'exit', figsize=(6, 3))

    #
    # Mouse
    #
    print('Plotting Mouse')
    mm_entry_deg = get_diff_expr_genes('MM', 'Cyte_vs_Gonia', 'Up')
    mm_exit_deg = get_diff_expr_genes('MM', 'Tid_vs_Cyte', 'Down')

    df_mm = pd.read_csv('../01-diff-gene-exp/data/DGE/MM/MM_RawCounts_AllSperm.csv', index_col=0)
    #
    dict_cols_gonia = {
        'Spermatogonia_Rep_1': 'Gonia-1',
        'Spermatogonia_Rep_2': 'Gonia-2',
        'Spermatogonia_Rep_3': 'Gonia-3',
    }
    dict_cols_cyte = {
        'Spermatocytes_Rep_1': 'Cyte-1',
        'Spermatocytes_Rep_2': 'Cyte-2',
        'Spermatocytes_Rep_3': 'Cyte-3',
    }
    dict_cols_tid = {
        'Spermatids_Rep_1': 'Tid-1',
        'Spermatids_Rep_2': 'Tid-2',
        'Spermatids_Rep_3': 'Tid-3',
    }
    dict_all_cols = {**dict_cols_gonia, **dict_cols_cyte, **dict_cols_tid}
    cols_gonia = list(dict_cols_gonia.values())
    cols_cyte = list(dict_cols_cyte.values())
    cols_tid = list(dict_cols_tid.values())
    df_mm = df_mm.rename(columns=dict_all_cols)
    #
    df_mm_entry = df_mm.loc[df_mm.index.isin(mm_entry_deg), cols_gonia + cols_cyte]
    df_mm_exit = df_mm.loc[df_mm.index.isin(mm_exit_deg), cols_cyte + cols_tid]
    #
    #df_mm_entry = pd.DataFrame(preprocessing.normalize(df_mm_entry, norm='l2'), columns=df_mm_entry.columns, index=df_mm_entry.index)
    #df_mm_exit = pd.DataFrame(preprocessing.normalize(df_mm_exit, norm='l2'), columns=df_mm_exit.columns, index=df_mm_exit.index)
    df_mm_entry = probnorm(df_mm_entry)
    df_mm_exit = probnorm(df_mm_exit)
    #
    plot(df_mm_entry, 'MM', 'entry', figsize=(3.8, 3))
    plot(df_mm_exit, 'MM', 'exit', figsize=(3.8, 3))

    #
    # Insect
    #
    print('Plotting Insect')
    dm_entry_deg = get_diff_expr_genes('DM', 'Cyte_vs_Gonia', 'Up')
    dm_exit_deg = get_diff_expr_genes('DM', 'Tid_vs_Cyte', 'Down')

    df_dm = pd.read_csv('../01-diff-gene-exp/data/DGE/DM/DM_RawCounts_AllSperm.csv', index_col=0)
    #
    dict_cols_apical = {
        'Apical_1': 'Apical-1',
        'Apical_2':' Apical-2',
    }
    dict_cols_middle = {
        'Middle_1': 'Middle-1',
        'Middle_2': 'Middle-2',
    }
    dict_cols_basal = {
        'Basal_1': 'Basal-1',
        'Basal_2': 'Basal-2',
    }
    dict_all_cols = {**dict_cols_apical, **dict_cols_middle, **dict_cols_basal}
    cols_apical = list(dict_cols_apical.values())
    cols_middle = list(dict_cols_middle.values())
    cols_basal = list(dict_cols_basal.values())
    df_dm = df_dm.rename(columns=dict_all_cols)
    #
    df_dm_entry = df_dm.loc[df_dm.index.isin(dm_entry_deg), cols_apical + cols_middle]
    df_dm_exit = df_dm.loc[df_dm.index.isin(dm_exit_deg), cols_middle + cols_basal]
    #
    #df_dm_entry = pd.DataFrame(preprocessing.normalize(df_dm_entry, norm='l2'), columns=df_dm_entry.columns, index=df_dm_entry.index)
    #df_dm_exit = pd.DataFrame(preprocessing.normalize(df_dm_exit, norm='l2'), columns=df_dm_exit.columns, index=df_dm_exit.index)
    df_dm_entry = probnorm(df_dm_entry)
    df_dm_exit = probnorm(df_dm_exit)
    #
    plot(df_dm_entry, 'DM', 'entry', figsize=(3, 3))
    plot(df_dm_exit, 'DM', 'exit', figsize=(3, 3))
