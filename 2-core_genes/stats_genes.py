# coding=utf-8
# Author: Rion B Correia
# Date: Aug 06, 2019
#
# Description: Dispalys Statistics about number of genes (including DGE) we found for each species
#
# Instructions:
#
import math
# import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.precision', 4)
from tabulate import tabulate


def df2md(df, y_index=False, *args, **kwargs):
    blob = tabulate(df, headers='keys', tablefmt='pipe', *args, **kwargs)
    if not y_index:
        return '\n'.join(['| {}'.format(row.split('|', 2)[-1]) for row in blob.split('\n')])
    return blob


if __name__ == '__main__':

    maxFDR = 0.05
    minLogFC = math.log2(2)

    #
    # HS
    #
    df_HS = pd.read_csv('results/HS-DE_genes.csv', index_col=0)

    df_HS_CyteGonia = df_HS.loc[(df_HS['Cyte_vs_Gonia']) == True, :]
    df_HS_TidCyte = df_HS.loc[(df_HS['Tid_vs_Cyte']) == True, :]

    # DE (Up/Down/Not)
    df_HS_DE_UpCyteGonia = df_HS_CyteGonia.loc[((df_HS_CyteGonia['FDR_CyteGonia'] <= maxFDR) & (df_HS_CyteGonia['logFC_CyteGonia'].abs() >= minLogFC) & (df_HS_CyteGonia['logFC_CyteGonia'] >= 0)), :]
    df_HS_DE_DownCyteGonia = df_HS_CyteGonia.loc[((df_HS_CyteGonia['FDR_CyteGonia'] <= maxFDR) & (df_HS_CyteGonia['logFC_CyteGonia'].abs() >= minLogFC) & (df_HS_CyteGonia['logFC_CyteGonia'] <= 0)), :]
    df_HS_DE_NotCyteGonia = df_HS_CyteGonia.loc[~df_HS_CyteGonia.index.isin(df_HS_DE_UpCyteGonia.index.tolist() + df_HS_DE_DownCyteGonia.index.tolist()), :]
    df_HS_DE_UpCyteTid = df_HS_TidCyte.loc[((df_HS_TidCyte['FDR_TidCyte'] <= maxFDR) & (df_HS_TidCyte['logFC_TidCyte'].abs() >= minLogFC) & (df_HS_TidCyte['logFC_TidCyte'] >= 0)), :]
    df_HS_DE_DownCyteTid = df_HS_TidCyte.loc[((df_HS_TidCyte['FDR_TidCyte'] <= maxFDR) & (df_HS_TidCyte['logFC_TidCyte'].abs() >= minLogFC) & (df_HS_TidCyte['logFC_TidCyte'] <= 0)), :]
    df_HS_DE_NotCyteTid = df_HS_TidCyte.loc[~df_HS_TidCyte.index.isin(df_HS_DE_UpCyteTid.index.tolist() + df_HS_DE_DownCyteTid.index.tolist()), :]

    #
    # HS (Mitosis)
    #
    df_HS_mit = pd.read_csv('results/HS-E_mitotic_genes.csv', index_col=0)
    df_HS_MitPre = df_HS_mit.loc[(df_HS_mit['Mit_vs_Pre']) == True, :]
    df_HS_PosMit = df_HS_mit.loc[(df_HS_mit['Pos_vs_Mit']) == True, :]

    #
    # MM
    #
    df_MM = pd.read_csv('results/MM-DE_genes.csv', index_col=0)

    df_MM_CyteGonia = df_MM.loc[(df_MM['Cyte_vs_Gonia']) == True, :]
    df_MM_TidCyte = df_MM.loc[(df_MM['Tid_vs_Cyte']) == True, :]

    # DE (Up/Down/Not)
    df_MM_DE_UpCyteGonia = df_MM_CyteGonia.loc[((df_MM_CyteGonia['FDR_CyteGonia'] <= maxFDR) & (df_MM_CyteGonia['logFC_CyteGonia'].abs() >= minLogFC) & (df_MM_CyteGonia['logFC_CyteGonia'] >= 0)), :]
    df_MM_DE_DownCyteGonia = df_MM_CyteGonia.loc[((df_MM_CyteGonia['FDR_CyteGonia'] <= maxFDR) & (df_MM_CyteGonia['logFC_CyteGonia'].abs() >= minLogFC) & (df_MM_CyteGonia['logFC_CyteGonia'] <= 0)), :]
    df_MM_DE_NotCyteGonia = df_MM_CyteGonia.loc[~df_MM_CyteGonia.index.isin(df_MM_DE_UpCyteGonia.index.tolist() + df_MM_DE_DownCyteGonia.index.tolist()), :]
    df_MM_DE_UpCyteTid = df_MM_TidCyte.loc[((df_MM_TidCyte['FDR_TidCyte'] <= maxFDR) & (df_MM_TidCyte['logFC_TidCyte'].abs() >= minLogFC) & (df_MM_TidCyte['logFC_TidCyte'] >= 0)), :]
    df_MM_DE_DownCyteTid = df_MM_TidCyte.loc[((df_MM_TidCyte['FDR_TidCyte'] <= maxFDR) & (df_MM_TidCyte['logFC_TidCyte'].abs() >= minLogFC) & (df_MM_TidCyte['logFC_TidCyte'] <= 0)), :]
    df_MM_DE_NotCyteTid = df_MM_TidCyte.loc[~df_MM_TidCyte.index.isin(df_MM_DE_UpCyteTid.index.tolist() + df_MM_DE_DownCyteTid.index.tolist()), :]

    #
    # DM
    #
    df_DM = pd.read_csv('results/DM-DE_genes.csv', index_col=0)

    df_DM_MiddleApical = df_DM.loc[(df_DM['Middle_vs_Apical']) == True, :]
    df_DM_BasalMiddle = df_DM.loc[(df_DM['Basal_vs_Middle']) == True, :]

    # DE (Up/Down/Not)
    df_DM_DE_UpMiddleApical = df_DM_MiddleApical.loc[((df_DM_MiddleApical['FDR_MiddleApical'] <= maxFDR) & (df_DM_MiddleApical['logFC_MiddleApical'].abs() >= minLogFC) & (df_DM_MiddleApical['logFC_MiddleApical'] >= 0)), :]
    df_DM_DE_DownMiddleApical = df_DM_MiddleApical.loc[((df_DM_MiddleApical['FDR_MiddleApical'] <= maxFDR) & (df_DM_MiddleApical['logFC_MiddleApical'].abs() >= minLogFC) & (df_DM_MiddleApical['logFC_MiddleApical'] <= 0)), :]
    df_DM_DE_NotMiddleApical = df_DM_MiddleApical.loc[~df_DM_MiddleApical.index.isin(df_DM_DE_UpMiddleApical.index.tolist() + df_DM_DE_DownMiddleApical.index.tolist()), :]
    df_DM_DE_UpMiddleBasal = df_DM_BasalMiddle.loc[((df_DM_BasalMiddle['FDR_BasalMiddle'] <= maxFDR) & (df_DM_BasalMiddle['logFC_BasalMiddle'].abs() >= minLogFC) & (df_DM_BasalMiddle['logFC_BasalMiddle'] >= 0)), :]
    df_DM_DE_DownMiddleBasal = df_DM_BasalMiddle.loc[((df_DM_BasalMiddle['FDR_BasalMiddle'] <= maxFDR) & (df_DM_BasalMiddle['logFC_BasalMiddle'].abs() >= minLogFC) & (df_DM_BasalMiddle['logFC_BasalMiddle'] <= 0)), :]
    df_DM_DE_NotMiddleBasal = df_DM_BasalMiddle.loc[~df_DM_BasalMiddle.index.isin(df_DM_DE_UpMiddleBasal.index.tolist() + df_DM_DE_DownMiddleBasal.index.tolist()), :]

    #
    # - number of genes
    #
    n_HS_g = df_HS.shape[0]
    n_MM_g = df_MM.shape[0]
    n_DM_g = df_DM.shape[0]
    n_HS_CyteGonia_g = df_HS_CyteGonia.shape[0]
    n_HS_TidCyte_g = df_HS_TidCyte.shape[0]
    n_HS_MitPre = df_HS_MitPre.shape[0]
    n_HS_PosMit = df_HS_PosMit.shape[0]
    n_MM_CyteGonia_g = df_MM_CyteGonia.shape[0]
    n_MM_TidCyte_g = df_MM_TidCyte.shape[0]
    n_DM_MiddleApical_g = df_DM_MiddleApical.shape[0]
    n_DM_BasalMiddle_g = df_DM_BasalMiddle.shape[0]

    n_HS_pcg = df_HS['biotype'].value_counts()['protein_coding']
    n_MM_pcg = df_MM['biotype'].value_counts()['protein_coding']
    n_DM_pcg = df_DM['biotype'].value_counts()['protein_coding']
    n_HS_CyteGonia_pcg = df_HS_CyteGonia['biotype'].value_counts()['protein_coding']
    n_HS_TidCyte_pcg = df_HS_TidCyte['biotype'].value_counts()['protein_coding']
    n_HS_MitPre_pcg = df_HS_MitPre['biotype'].value_counts()['protein_coding']
    n_HS_PosMit_pcg = df_HS_PosMit['biotype'].value_counts()['protein_coding']
    n_MM_CyteGonia_pcg = df_MM_CyteGonia['biotype'].value_counts()['protein_coding']
    n_MM_TidCyte_pcg = df_MM_TidCyte['biotype'].value_counts()['protein_coding']
    n_DM_MiddleApical_pcg = df_DM_MiddleApical['biotype'].value_counts()['protein_coding']
    n_DM_BasalMiddle_pcg = df_DM_BasalMiddle['biotype'].value_counts()['protein_coding']

    print('# Number of genes\n')

    df_stat = pd.DataFrame.from_records([
        ('HS', 'Both', n_HS_g, n_HS_pcg),
        ('HS', 'Cyte vs Gonia', n_HS_CyteGonia_g, n_HS_CyteGonia_pcg),
        ('HS', 'Tid vs Cyte', n_HS_TidCyte_g, n_HS_TidCyte_pcg),
        ('HS (mitosis)', 'Mitosis vs Pre-Mitosis', n_HS_MitPre, n_HS_MitPre_pcg),
        ('HS (mitosis)', 'Pos-Mitosis vs Mitosis', n_HS_PosMit, n_HS_PosMit_pcg),
        ('MM', 'Both', n_MM_g, n_MM_pcg),
        ('MM', 'Cyte vs Gonia', n_MM_CyteGonia_g, n_MM_CyteGonia_pcg),
        ('MM', 'Tid vs Cyte', n_MM_TidCyte_g, n_MM_TidCyte_pcg),
        ('DM', 'Both', n_DM_g, n_DM_pcg),
        ('DM', 'Middle vs Apical', n_DM_MiddleApical_g, n_DM_MiddleApical_pcg),
        ('DM', 'Basal vs Middle', n_DM_BasalMiddle_g, n_DM_BasalMiddle_pcg),
    ], columns=['Species', 'Cell', 'Genes', 'Prot. Coding'])
    df_stat['%'] = df_stat['Prot. Coding'] / df_stat['Genes']
    print(df2md(df_stat, floatfmt='.4f'))
    print('\n')

    #
    # - number of genes differently expressed
    #
    print("# Number of genes differently expressed\n")

    # HS
    n_HS_UpCyteGonia_g = df_HS_DE_UpCyteGonia.shape[0]
    n_HS_DownCyteGonia_g = df_HS_DE_DownCyteGonia.shape[0]
    n_HS_NotCyteGonia_g = df_HS_DE_NotCyteGonia.shape[0]
    n_HS_UpCyteTid_g = df_HS_DE_UpCyteTid.shape[0]
    n_HS_DownCyteTid_g = df_HS_DE_DownCyteTid.shape[0]
    n_HS_NotCyteTid_g = df_HS_DE_NotCyteTid.shape[0]

    n_HS_UpCyteGonia_pcg = df_HS_DE_UpCyteGonia['biotype'].value_counts()['protein_coding']
    n_HS_DownCyteGonia_pcg = df_HS_DE_DownCyteGonia['biotype'].value_counts()['protein_coding']
    n_HS_NotCyteGonia_pcg = df_HS_DE_NotCyteGonia['biotype'].value_counts()['protein_coding']
    n_HS_UpCyteTid_pcg = df_HS_DE_UpCyteTid['biotype'].value_counts()['protein_coding']
    n_HS_DownCyteTid_pcg = df_HS_DE_DownCyteTid['biotype'].value_counts()['protein_coding']
    n_HS_NotCyteTid_pcg = df_HS_DE_NotCyteTid['biotype'].value_counts()['protein_coding']

    # MM
    n_MM_UpCyteGonia_g = df_MM_DE_UpCyteGonia.shape[0]
    n_MM_DownCyteGonia_g = df_MM_DE_DownCyteGonia.shape[0]
    n_MM_NotCyteGonia_g = df_MM_DE_NotCyteGonia.shape[0]
    n_MM_UpCyteTid_g = df_MM_DE_UpCyteTid.shape[0]
    n_MM_DownCyteTid_g = df_MM_DE_DownCyteTid.shape[0]
    n_MM_NotCyteTid_g = df_MM_DE_NotCyteTid.shape[0]

    n_MM_UpCyteGonia_pcg = df_MM_DE_UpCyteGonia['biotype'].value_counts()['protein_coding']
    n_MM_DownCyteGonia_pcg = df_MM_DE_DownCyteGonia['biotype'].value_counts()['protein_coding']
    n_MM_NotCyteGonia_pcg = df_MM_DE_NotCyteGonia['biotype'].value_counts()['protein_coding']
    n_MM_UpCyteTid_pcg = df_MM_DE_UpCyteTid['biotype'].value_counts()['protein_coding']
    n_MM_DownCyteTid_pcg = df_MM_DE_DownCyteTid['biotype'].value_counts()['protein_coding']
    n_MM_NotCyteTid_pcg = df_MM_DE_NotCyteTid['biotype'].value_counts()['protein_coding']

    # DM
    n_DM_UpMiddleApical_g = df_DM_DE_UpMiddleApical.shape[0]
    n_DM_DownMiddleApical_g = df_DM_DE_DownMiddleApical.shape[0]
    n_DM_NotMiddleApical_g = df_DM_DE_NotMiddleApical.shape[0]
    n_DM_UpMiddleBasal_g = df_DM_DE_UpMiddleBasal.shape[0]
    n_DM_DownMiddleBasal_g = df_DM_DE_DownMiddleBasal.shape[0]
    n_DM_NotMiddleBasal_g = df_DM_DE_NotMiddleBasal.shape[0]

    n_DM_UpMiddleApical_pcg = df_DM_DE_UpMiddleApical['biotype'].value_counts()['protein_coding']
    n_DM_DownMiddleApical_pcg = df_DM_DE_DownMiddleApical['biotype'].value_counts()['protein_coding']
    n_DM_NotMiddleApical_pcg = df_DM_DE_NotMiddleApical['biotype'].value_counts()['protein_coding']
    n_DM_UpMiddleBasal_pcg = df_DM_DE_UpMiddleBasal['biotype'].value_counts()['protein_coding']
    n_DM_DownMiddleBasal_pcg = df_DM_DE_DownMiddleBasal['biotype'].value_counts()['protein_coding']
    n_DM_NotMiddleBasal_pcg = df_DM_DE_NotMiddleBasal['biotype'].value_counts()['protein_coding']

    columns = ['Specie', 'Cell', 'Reg.', 'FDR', 'Genes']

    # HS CyteGonia
    df_HS_CyteGonia_stat = pd.DataFrame.from_records([
        ('HS', 'Cyte vs Gonia', 'Up', maxFDR, n_HS_UpCyteGonia_g),
        ('HS', 'Cyte vs Gonia', 'Not', maxFDR, n_HS_NotCyteGonia_g),
        ('HS', 'Cyte vs Gonia', 'Down', maxFDR, n_HS_DownCyteGonia_g)
    ], columns=columns)
    df_HS_CyteGonia_stat['%(G)'] = df_HS_CyteGonia_stat['Genes'] / df_HS_CyteGonia_stat['Genes'].sum()
    df_HS_CyteGonia_stat['Prot. Coding'] = [n_HS_UpCyteGonia_pcg, n_HS_NotCyteGonia_pcg, n_HS_DownCyteGonia_pcg]
    df_HS_CyteGonia_stat['%(PC)'] = df_HS_CyteGonia_stat['Prot. Coding'] / df_HS_CyteGonia_stat['Genes']

    # HS TidCyte
    df_HS_TidCyte_stat = pd.DataFrame.from_records([
        ('HS', 'Tid vs Cyte', 'Up', maxFDR, n_HS_UpCyteTid_g),
        ('HS', 'Tid vs Cyte', 'Not', maxFDR, n_HS_NotCyteTid_g),
        ('HS', 'Tid vs Cyte', 'Down', maxFDR, n_HS_DownCyteTid_g)
    ], columns=columns)
    df_HS_TidCyte_stat['%(G)'] = df_HS_TidCyte_stat['Genes'] / df_HS_TidCyte_stat['Genes'].sum()
    df_HS_TidCyte_stat['Prot. Coding'] = [n_HS_UpCyteTid_pcg, n_HS_NotCyteTid_pcg, n_HS_DownCyteTid_pcg]
    df_HS_TidCyte_stat['%(PC)'] = df_HS_TidCyte_stat['Prot. Coding'] / df_HS_TidCyte_stat['Genes']

    # MM CyteGonia
    df_MM_CyteGonia_stat = pd.DataFrame.from_records([
        ('MM', 'Cyte vs Gonia', 'Up', maxFDR, n_MM_UpCyteGonia_g),
        ('MM', 'Cyte vs Gonia', 'Not', maxFDR, n_MM_NotCyteGonia_g),
        ('MM', 'Cyte vs Gonia', 'Down', maxFDR, n_MM_DownCyteGonia_g)
    ], columns=columns)
    df_MM_CyteGonia_stat['%(G)'] = df_MM_CyteGonia_stat['Genes'] / df_MM_CyteGonia_stat['Genes'].sum()
    df_MM_CyteGonia_stat['Prot. Coding'] = [n_MM_UpCyteGonia_pcg, n_MM_NotCyteGonia_pcg, n_MM_DownCyteGonia_pcg]
    df_MM_CyteGonia_stat['%(PC)'] = df_MM_CyteGonia_stat['Prot. Coding'] / df_MM_CyteGonia_stat['Genes']

    # MM TidCyte
    df_MM_TidCyte_stat = pd.DataFrame.from_records([
        ('MM', 'Tid vs Cyte', 'Up', maxFDR, n_MM_UpCyteTid_g),
        ('MM', 'Tid vs Cyte', 'Not', maxFDR, n_MM_NotCyteTid_g),
        ('MM', 'Tid vs Cyte', 'Down', maxFDR, n_MM_DownCyteTid_g)
    ], columns=columns)
    df_MM_TidCyte_stat['%(G)'] = df_MM_TidCyte_stat['Genes'] / df_MM_TidCyte_stat['Genes'].sum()
    df_MM_TidCyte_stat['Prot. Coding'] = [n_MM_UpCyteTid_pcg, n_MM_NotCyteTid_pcg, n_MM_DownCyteTid_pcg]
    df_MM_TidCyte_stat['%(PC)'] = df_MM_TidCyte_stat['Prot. Coding'] / df_MM_TidCyte_stat['Genes']

    # DM MiddleApical
    df_DM_MiddleApical_stat = pd.DataFrame.from_records([
        ('DM', 'Middle vs Apical', 'Up', maxFDR, n_DM_UpMiddleApical_g),
        ('DM', 'Middle vs Apical', 'Not', maxFDR, n_DM_NotMiddleApical_g),
        ('DM', 'Middle vs Apical', 'Down', maxFDR, n_DM_DownMiddleApical_g)
    ], columns=columns)
    df_DM_MiddleApical_stat['%(G)'] = df_DM_MiddleApical_stat['Genes'] / df_DM_MiddleApical_stat['Genes'].sum()
    df_DM_MiddleApical_stat['Prot. Coding'] = [n_DM_UpMiddleApical_pcg, n_DM_NotMiddleApical_pcg, n_DM_DownMiddleApical_pcg]
    df_DM_MiddleApical_stat['%(PC)'] = df_DM_MiddleApical_stat['Prot. Coding'] / df_DM_MiddleApical_stat['Genes']

    # DM BasalMiddle
    df_DM_BasalMiddle_stat = pd.DataFrame.from_records([
        ('DM', 'Basal vs Middle', 'Up', maxFDR, n_DM_UpMiddleBasal_g),
        ('DM', 'Basal vs Middle', 'Not', maxFDR, n_DM_NotMiddleBasal_g),
        ('DM', 'Basal vs Middle', 'Down', maxFDR, n_DM_DownMiddleBasal_g)
    ], columns=columns)
    df_DM_BasalMiddle_stat['%(G)'] = df_DM_BasalMiddle_stat['Genes'] / df_DM_BasalMiddle_stat['Genes'].sum()
    df_DM_BasalMiddle_stat['Prot. Coding'] = [n_DM_UpMiddleBasal_pcg, n_DM_NotMiddleBasal_pcg, n_DM_DownMiddleBasal_pcg]
    df_DM_BasalMiddle_stat['%(PC)'] = df_DM_BasalMiddle_stat['Prot. Coding'] / df_DM_BasalMiddle_stat['Genes']

    floatfmt = ['', '', '', '', '.2f', '', '.4f', '', '.4f']
    print(df2md(df_HS_CyteGonia_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_HS_TidCyte_stat, floatfmt=floatfmt))
    print('')
    #
    print(df2md(df_MM_CyteGonia_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_MM_TidCyte_stat, floatfmt=floatfmt))
    print('')
    #
    print(df2md(df_DM_MiddleApical_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_DM_BasalMiddle_stat, floatfmt=floatfmt))
    print('')
