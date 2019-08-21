# coding=utf-8
# Author: Rion B Correia
# Date: Aug 06, 2019
#
# Description: Dispalys Statistics about number of genes (including DGE) we found for each species
#
# Instructions:
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from tabulate import tabulate

def df2md(df, y_index=False):
    blob = tabulate(df, headers='keys', tablefmt='pipe')
    if not y_index:
        return '\n'.join(['| {}'.format(row.split('|', 2)[-1]) for row in blob.split('\n')])
    return blob

if __name__ == '__main__':

    #
    # - number of genes
    #
    print('# Number of genes\n')

    # HS     
    df_HS_DGE_CyteGonia = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df_HS_DGE_CyteTid = pd.read_csv('../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Tid.csv', index_col=0)
    # Fix Gene ID
    df_HS_DGE_CyteGonia.index = df_HS_DGE_CyteGonia.index.map(lambda x: x.split('.')[0])
    df_HS_DGE_CyteTid.index = df_HS_DGE_CyteTid.index.map(lambda x: x.split('.')[0])

    # MM
    df_MM_DGE_CyteGonia = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df_MM_DGE_CyteTid = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Tid.csv', index_col=0)

    # DM
    df_DM_DGE_Apical = pd.read_csv('../1-diff-gene-exp/results/DM/DM_DGE_ApicalTestis.csv', index_col=0)
    df_DM_DGE_MidTestis = pd.read_csv('../1-diff-gene-exp/results/DM/DM_DGE_MidTestis.csv', index_col=0)

    print(df2md(
        pd.DataFrame.from_records([
            ('HS', 'Cyte vs Gonia', np.unique(df_HS_DGE_CyteGonia.index).shape[0]),
            ('HS', 'Cyte vs Tid', np.unique(df_HS_DGE_CyteTid.index).shape[0]),
            ('MM', 'Cyte vs Gonia', np.unique(df_MM_DGE_CyteGonia.index).shape[0]),
            ('MM', 'Cyte vs Tid', np.unique(df_MM_DGE_CyteTid.index).shape[0]),
            ('DM', 'Apical Testis', np.unique(df_DM_DGE_Apical.index).shape[0]),
            ('DM', 'Mid Testis', np.unique(df_DM_DGE_MidTestis.index).shape[0]),
        ], columns=['Species', 'Cell', '# Genes'])
    ))

    print('\n')
    del df_HS_DGE_CyteGonia, df_HS_DGE_CyteTid,
    del df_MM_DGE_CyteGonia, df_MM_DGE_CyteTid,
    del df_DM_DGE_Apical, df_DM_DGE_MidTestis
    #
    # - number of genes differently expressed
    #
    print("# Number of genes differently expressed\n")
    
    # HS     
    df_HS_DE_UpCyteGonia01 = pd.read_csv('../1-diff-gene-exp/results/HS/DE/HS_DE_UpCyte_vs_Gonia-FDR_0p01.csv', index_col=0)
    df_HS_DE_UpCyteGonia05 = pd.read_csv('../1-diff-gene-exp/results/HS/DE/HS_DE_UpCyte_vs_Gonia-FDR_0p05.csv', index_col=0)
    df_HS_DE_DownCyteTid01 = pd.read_csv('../1-diff-gene-exp/results/HS/DE/HS_DE_DownCyte_vs_Tid-FDR_0p01.csv', index_col=0)
    df_HS_DE_DownCyteTid05 = pd.read_csv('../1-diff-gene-exp/results/HS/DE/HS_DE_DownCyte_vs_Tid-FDR_0p05.csv', index_col=0)

    # MM
    df_MM_DE_UpCyteGonia01 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_UpCyte_vs_Gonia-FDR_0p01.csv', index_col=0)
    df_MM_DE_UpCyteGonia05 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_UpCyte_vs_Gonia-FDR_0p05.csv', index_col=0)
    df_MM_DE_DownCyteTid01 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_DownCyte_vs_Tid-FDR_0p01.csv', index_col=0)
    df_MM_DE_DownCyteTid05 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_DownCyte_vs_Tid-FDR_0p05.csv', index_col=0)

    # DM
    df_DM_DE_UpApical = pd.read_csv('../1-diff-gene-exp/results/DM/DE/DM_DE_UpApicalTestis.csv', index_col=0)
    df_DM_DE_DownMid = pd.read_csv('../1-diff-gene-exp/results/DM/DE/DM_DE_DownMidTestis.csv', index_col=0)

    print(df2md(
        pd.DataFrame.from_records([
            ('HS', 'Up', 'Cyte vs Gonia', '<=0.01', np.unique(df_HS_DE_UpCyteGonia01.index).shape[0]),
            ('HS', 'Up', 'Cyte vs Gonia', '<=0.05', np.unique(df_HS_DE_UpCyteGonia05.index).shape[0]),
            ('HS', 'Down', 'Cyte vs Tid', '<=0.01', np.unique(df_HS_DE_DownCyteTid01.index).shape[0]),
            ('HS', 'Down', 'Cyte vs Tid', '<=0.05', np.unique(df_HS_DE_DownCyteTid05.index).shape[0]),
            #
            ('MM', 'Up', 'Cyte vs Gonia', '<=0.01', np.unique(df_MM_DE_UpCyteGonia01.index).shape[0]),
            ('MM', 'Up', 'Cyte vs Gonia', '<=0.05', np.unique(df_MM_DE_UpCyteGonia05.index).shape[0]),
            ('MM', 'Down', 'Cyte vs Tid', '<=0.01', np.unique(df_MM_DE_DownCyteTid01.index).shape[0]),
            ('MM', 'Down', 'Cyte vs Tid', '<=0.05', np.unique(df_MM_DE_DownCyteTid05.index).shape[0]),
            #
            ('DM', 'Up', 'Apical Testis', 'n/a', np.unique(df_DM_DE_UpApical.index).shape[0]),
            ('DM', 'Down', 'Mid Testis', 'n/a', np.unique(df_DM_DE_DownMid.index).shape[0]),
        ], columns=['Species', 'Regulation', 'Cell', 'FDR', '# Genes'])
    ))
    print('\n')

    del df_HS_DE_UpCyteGonia01, df_HS_DE_UpCyteGonia05, df_HS_DE_DownCyteTid01, df_HS_DE_DownCyteTid05
    del df_MM_DE_UpCyteGonia01, df_MM_DE_UpCyteGonia05, df_MM_DE_DownCyteTid01, df_MM_DE_DownCyteTid05
    del df_DM_DE_UpApical, df_DM_DE_DownMid

    #
    # - number of protein-coding genes
    #

    #
    # - number of protein-coding genes differently expressed
    #
    print('# Number of protein-coding genes differently expressed\n')

    # HS
    df_HS_genes01 = pd.read_csv('results/HS/genes_HS-FDR_0p01.csv', index_col=0)
    df_HS_genes05 = pd.read_csv('results/HS/genes_HS-FDR_0p05.csv', index_col=0)

    print('## HS with FDR=0.01\n')
    print(df2md(df_HS_genes01['biotype'].value_counts().to_frame(), y_index=True))
    print('\n## HS with: FDR=0.05\n')
    print(df2md(df_HS_genes05['biotype'].value_counts().to_frame(), y_index=True))

    # MM
    df_MM_genes01 = pd.read_csv('results/MM/genes_MM-FDR_0p01.csv', index_col=0)
    df_MM_genes05 = pd.read_csv('results/MM/genes_MM-FDR_0p05.csv', index_col=0)

    print('\n## MM with FDR=0.01\n')
    print(df2md(df_MM_genes01['biotype'].value_counts().to_frame(), y_index=True))
    print('\n## MM with FDR=0.05\n')
    print(df2md(df_MM_genes05['biotype'].value_counts().to_frame(), y_index=True))

    # DS
    df_DS_genes = pd.read_csv('results/DM/genes_DM.csv', index_col=0)

    print('\n## DS\n')
    print(df2md(df_DS_genes['biotype'].value_counts().to_frame(), y_index=True))

    #
    # - number of genes found conserved in other species (only for the gene populations that were compared across species)
    #
