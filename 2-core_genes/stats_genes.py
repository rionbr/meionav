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

def dataframe2markdown(df, y_index=False):
    blob = tabulate(df, headers='keys', tablefmt='pipe')
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

    print('- HS Cyte vs Gonia: {:,d}'.format(np.unique(df_HS_DGE_CyteGonia.index).shape[0]))
    print('- HS Cyte vs Tid: {:,d}\n'.format(np.unique(df_HS_DGE_CyteTid.index).shape[0]))

    # MM
    df_MM_DGE_CyteGonia = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df_MM_DGE_CyteTid = pd.read_csv('../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Tid.csv', index_col=0)

    print('- MM Cyte vs Gonia: {:,d}'.format(np.unique(df_MM_DGE_CyteGonia.index).shape[0]))
    print('- MM Cyte vs Tid: {:,d}\n'.format(np.unique(df_MM_DGE_CyteTid.index).shape[0]))
    
    # DM
    df_DM_DGE_Apical = pd.read_csv('../1-diff-gene-exp/results/DM/DM_DGE_ApicalTestis.csv', index_col=0)
    df_DM_DGE_MidTestis = pd.read_csv('../1-diff-gene-exp/results/DM/DM_DGE_MidTestis.csv', index_col=0)

    print('- DM Apical: {:,d}'.format(np.unique(df_DM_DGE_Apical.index).shape[0]))
    print('- DM MidTestis: {:,d}\n'.format(np.unique(df_DM_DGE_MidTestis.index).shape[0]))

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

    print('- MM (Up)Cyte vs Gonia: FDR=0.01: {:,d}; FDR=0.05: {:,d}'.format(
        np.unique(df_HS_DE_UpCyteGonia01.index).shape[0],
        np.unique(df_HS_DE_UpCyteGonia05.index).shape[0]))
    print('- MM (Down)Cyte vs Tid: FDR=0.01: {:,d}; FDR=0.05: {:,d}\n'.format(
        np.unique(df_HS_DE_DownCyteTid01.index).shape[0],
        np.unique(df_HS_DE_DownCyteTid05.index).shape[0]))

    # MM
    df_MM_DE_UpCyteGonia01 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_UpCyte_vs_Gonia-FDR_0p01.csv', index_col=0)
    df_MM_DE_UpCyteGonia05 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_UpCyte_vs_Gonia-FDR_0p05.csv', index_col=0)
    df_MM_DE_DownCyteTid01 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_DownCyte_vs_Tid-FDR_0p01.csv', index_col=0)
    df_MM_DE_DownCyteTid05 = pd.read_csv('../1-diff-gene-exp/results/MM/DE/MM_DE_DownCyte_vs_Tid-FDR_0p05.csv', index_col=0)

    print('- MM (Up)Cyte vs Gonia: FDR=0.01: {:,d}; FDR=0.05: {:,d}'.format(
        np.unique(df_MM_DE_UpCyteGonia01.index).shape[0],
        np.unique(df_MM_DE_UpCyteGonia05.index).shape[0]))
    print('- MM (Down)Cyte vs Tid: FDR=0.01: {:,d}; FDR=0.05: {:,d}\n'.format(
        np.unique(df_MM_DE_DownCyteTid01.index).shape[0],
        np.unique(df_MM_DE_DownCyteTid05.index).shape[0]))

    # DM
    df_DM_DE_UpApical = pd.read_csv('../1-diff-gene-exp/results/DM/DE/DM_DE_UpApicalTestis.csv', index_col=0)
    df_DM_DE_DownMid = pd.read_csv('../1-diff-gene-exp/results/DM/DE/DM_DE_DownMidTestis.csv', index_col=0)

    print('- DM (Up)Apical Testis: {:,d}'.format(np.unique(df_DM_DE_UpApical.index).shape[0]))
    print('- DM (Down)MidTestis: {:,d}\n'.format(np.unique(df_DM_DE_DownMid.index).shape[0]))
    
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

    print('## HS with FDR=0.01:\n')
    print(dataframe2markdown(df_HS_genes01['biotype'].value_counts().to_frame()))
    print('\n## HS with: FDR=0.05:\n')
    print(dataframe2markdown(df_HS_genes05['biotype'].value_counts().to_frame()))

    # MM
    df_MM_genes01 = pd.read_csv('results/MM/genes_MM-FDR_0p01.csv', index_col=0)
    df_MM_genes05 = pd.read_csv('results/MM/genes_MM-FDR_0p05.csv', index_col=0)

    print('\n## MM with FDR=0.01:\n')
    print(dataframe2markdown(df_MM_genes01['biotype'].value_counts().to_frame()))
    print('\n## MM with FDR=0.05:\n')
    print(dataframe2markdown(df_MM_genes05['biotype'].value_counts().to_frame()))

    # DS
    df_DS_genes = pd.read_csv('results/DM/genes_DM.csv', index_col=0)

    print('\n## DS:\n')
    print(dataframe2markdown(df_DS_genes['biotype'].value_counts().to_frame()))

    #
    # - number of genes found conserved in other species (only for the gene populations that were compared across species)
    #
