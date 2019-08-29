# coding=utf-8
# Author: Rion B Correia
# Date: Aug 06, 2019
#
# Description: Dispalys Statistics about number of genes (including DGE) we found for each species
#
# Instructions:
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.precision', 4)
from tabulate import tabulate
from pybiomart import Dataset


def df2md(df, y_index=False, *args, **kwargs):
    blob = tabulate(df, headers='keys', tablefmt='pipe', *args, **kwargs)
    if not y_index:
        return '\n'.join(['| {}'.format(row.split('|', 2)[-1]) for row in blob.split('\n')])
    return blob

if __name__ == '__main__':

    maxFDR = 0.05
    minLogFC = math.log2(2)

    # Query bioMart for Gene Name/Description    
    ds_HS = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    df_HS_G = ds_HS.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')
    ds_MM = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
    df_MM_G = ds_MM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')
    ds_DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    df_DM_G = ds_DM.query(attributes=['affy_drosophila_2', 'ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID')
    #
    # HS 
    #
    df_HS_DGE_CyteGonia = pd.read_csv('results/HS/HS-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df_HS_DGE_CyteTid = pd.read_csv('results/HS/HS-DGE_Cyte_vs_Tid.csv', index_col=0)
    # Fix Gene ID
    df_HS_DGE_CyteGonia.index = df_HS_DGE_CyteGonia.index.map(lambda x: x.split('.')[0])
    df_HS_DGE_CyteTid.index = df_HS_DGE_CyteTid.index.map(lambda x: x.split('.')[0])
    # Add Additional Info
    df_HS_DGE_CyteGonia['biotype'] = df_HS_G['Gene type']
    df_HS_DGE_CyteTid['biotype'] = df_HS_G['Gene type']
    # DE (Up/Down/Not)
    df_HS_DE_UpCyteGonia = df_HS_DGE_CyteGonia.loc[((df_HS_DGE_CyteGonia['FDR'] <= maxFDR) & (df_HS_DGE_CyteGonia['logFC'].abs() >= minLogFC) & (df_HS_DGE_CyteGonia['logFC'] >= 0)), :]
    df_HS_DE_DownCyteGonia = df_HS_DGE_CyteGonia.loc[((df_HS_DGE_CyteGonia['FDR'] <= maxFDR) & (df_HS_DGE_CyteGonia['logFC'].abs() >= minLogFC) & (df_HS_DGE_CyteGonia['logFC'] <= 0)), :]
    df_HS_DE_NotCyteGonia = df_HS_DGE_CyteGonia.loc[~df_HS_DGE_CyteGonia.index.isin(df_HS_DE_UpCyteGonia.index.tolist() + df_HS_DE_DownCyteGonia.index.tolist()), :]
    df_HS_DE_UpCyteTid = df_HS_DGE_CyteTid.loc[((df_HS_DGE_CyteTid['FDR'] <= maxFDR) & (df_HS_DGE_CyteTid['logFC'].abs() >= minLogFC) & (df_HS_DGE_CyteTid['logFC'] >= 0)), :]
    df_HS_DE_DownCyteTid = df_HS_DGE_CyteTid.loc[((df_HS_DGE_CyteTid['FDR'] <= maxFDR) & (df_HS_DGE_CyteTid['logFC'].abs() >= minLogFC) & (df_HS_DGE_CyteTid['logFC'] <= 0)), :]
    df_HS_DE_NotCyteTid = df_HS_DGE_CyteTid.loc[~df_HS_DGE_CyteTid.index.isin(df_HS_DE_UpCyteTid.index.tolist() + df_HS_DE_DownCyteTid.index.tolist()), :]

    #
    # MM
    #
    df_MM_DGE_CyteGonia = pd.read_csv('results/MM/MM-DGE_Cyte_vs_Gonia.csv', index_col=0)
    df_MM_DGE_CyteTid = pd.read_csv('results/MM/MM-DGE_Cyte_vs_Tid.csv', index_col=0)
    # Add Additional Info
    df_MM_DGE_CyteGonia['biotype'] = df_MM_G['Gene type']
    df_MM_DGE_CyteTid['biotype'] = df_MM_G['Gene type']
    # DE (Up/Down/Not)
    df_MM_DE_UpCyteGonia = df_MM_DGE_CyteGonia.loc[((df_MM_DGE_CyteGonia['FDR'] <= maxFDR) & (df_MM_DGE_CyteGonia['logFC'].abs() >= minLogFC) & (df_MM_DGE_CyteGonia['logFC'] >= 0)), :]
    df_MM_DE_DownCyteGonia = df_MM_DGE_CyteGonia.loc[((df_MM_DGE_CyteGonia['FDR'] <= maxFDR) & (df_MM_DGE_CyteGonia['logFC'].abs() >= minLogFC) & (df_MM_DGE_CyteGonia['logFC'] <= 0)), :]
    df_MM_DE_NotCyteGonia = df_MM_DGE_CyteGonia.loc[~df_MM_DGE_CyteGonia.index.isin(df_MM_DE_UpCyteGonia.index.tolist() + df_MM_DE_DownCyteGonia.index.tolist()), :]
    df_MM_DE_UpCyteTid = df_MM_DGE_CyteTid.loc[((df_MM_DGE_CyteTid['FDR'] <= maxFDR) & (df_MM_DGE_CyteTid['logFC'].abs() >= minLogFC) & (df_MM_DGE_CyteTid['logFC'] >= 0)), :]
    df_MM_DE_DownCyteTid = df_MM_DGE_CyteTid.loc[((df_MM_DGE_CyteTid['FDR'] <= maxFDR) & (df_MM_DGE_CyteTid['logFC'].abs() >= minLogFC) & (df_MM_DGE_CyteTid['logFC'] <= 0)), :]
    df_MM_DE_NotCyteTid = df_MM_DGE_CyteTid.loc[~df_MM_DGE_CyteTid.index.isin(df_MM_DE_UpCyteTid.index.tolist() + df_MM_DE_DownCyteTid.index.tolist()), :]
    
    #
    # DM
    #
    df_DM_DGE_MiddleApical = pd.read_csv('results/DM/DM-DGE_Middle_vs_Apical.csv', index_col=0)
    df_DM_DGE_MiddleBasal = pd.read_csv('results/DM/DM-DGE_Middle_vs_Basal.csv', index_col=0)
    
    df_DM_DGE_MiddleApical = pd.merge(df_DM_DGE_MiddleApical, df_DM_G, left_index=True, right_index=True, how='left')
    df_DM_DGE_MiddleBasal = pd.merge(df_DM_DGE_MiddleBasal, df_DM_G, left_index=True, right_index=True, how='left')

    # Note we rename column 'adj.P.Val' for 'FDR' for simplicity.
    df_DM_DGE_MiddleApical.rename(columns={'adj.P.Val': 'FDR', 'Gene type': 'biotype'}, inplace=True)
    df_DM_DGE_MiddleBasal.rename(columns={'adj.P.Val': 'FDR', 'Gene type': 'biotype'}, inplace=True)

    # DE (Up/Down/Not)
    df_DM_DE_UpMiddleApical = df_DM_DGE_MiddleApical.loc[((df_DM_DGE_MiddleApical['FDR'] <= maxFDR) & (df_DM_DGE_MiddleApical['logFC'].abs() >= minLogFC) & (df_DM_DGE_MiddleApical['logFC'] >= 0)), :]
    df_DM_DE_DownMiddleApical = df_DM_DGE_MiddleApical.loc[((df_DM_DGE_MiddleApical['FDR'] <= maxFDR) & (df_DM_DGE_MiddleApical['logFC'].abs() >= minLogFC) & (df_DM_DGE_MiddleApical['logFC'] <= 0)), :]
    df_DM_DE_NotMiddleApical = df_DM_DGE_MiddleApical.loc[~df_DM_DGE_MiddleApical.index.isin(df_DM_DE_UpMiddleApical.index.tolist() + df_DM_DE_DownMiddleApical.index.tolist()), :]
    df_DM_DE_UpMiddleBasal = df_DM_DGE_MiddleBasal.loc[((df_DM_DGE_MiddleBasal['FDR'] <= maxFDR) & (df_DM_DGE_MiddleBasal['logFC'].abs() >= minLogFC) & (df_DM_DGE_MiddleBasal['logFC'] >= 0)), :]
    df_DM_DE_DownMiddleBasal = df_DM_DGE_MiddleBasal.loc[((df_DM_DGE_MiddleBasal['FDR'] <= maxFDR) & (df_DM_DGE_MiddleBasal['logFC'].abs() >= minLogFC) & (df_DM_DGE_MiddleBasal['logFC'] <= 0)), :]
    df_DM_DE_NotMiddleBasal = df_DM_DGE_MiddleBasal.loc[~df_DM_DGE_MiddleBasal.index.isin(df_DM_DE_UpMiddleBasal.index.tolist() + df_DM_DE_DownMiddleBasal.index.tolist()), :]

    #
    # - number of genes
    #
    n_HS_CyteGonia_g = df_HS_DGE_CyteGonia.shape[0]
    n_HS_CyteTid_g = df_HS_DGE_CyteTid.shape[0]
    n_MM_CyteGonia_g = df_MM_DGE_CyteGonia.shape[0]
    n_MM_CyteTid_g = df_MM_DGE_CyteTid.shape[0]
    n_DM_MiddleApical_g = df_DM_DGE_MiddleApical.shape[0]
    n_DM_MiddleBasal_g = df_DM_DGE_MiddleBasal.shape[0]

    n_HS_CyteGonia_pcg = df_HS_DGE_CyteGonia['biotype'].value_counts()['protein_coding']
    n_HS_CyteTid_pcg = df_HS_DGE_CyteTid['biotype'].value_counts()['protein_coding']
    n_MM_CyteGonia_pcg = df_MM_DGE_CyteGonia['biotype'].value_counts()['protein_coding']
    n_MM_CyteTid_pcg = df_MM_DGE_CyteTid['biotype'].value_counts()['protein_coding']
    n_DM_MiddleApical_pcg = df_DM_DGE_MiddleApical['biotype'].value_counts()['protein_coding']
    n_DM_MiddleBasal_pcg = df_DM_DGE_MiddleBasal['biotype'].value_counts()['protein_coding']

    print('# Number of genes\n')

    df_stat = pd.DataFrame.from_records([
        ('HS', 'Cyte vs Gonia', n_HS_CyteGonia_g, n_HS_CyteGonia_pcg),
        ('HS', 'Cyte vs Tid', n_HS_CyteTid_g, n_HS_CyteTid_pcg),
        ('MM', 'Cyte vs Gonia', n_MM_CyteGonia_g, n_MM_CyteGonia_pcg),
        ('MM', 'Cyte vs Tid', n_MM_CyteTid_g, n_MM_CyteTid_pcg),
        ('DM', 'Middle vs Apical', n_DM_MiddleApical_g, n_DM_MiddleApical_pcg),
        ('DM', 'Middle vs Basal', n_DM_MiddleBasal_g, n_DM_MiddleBasal_pcg),
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
    
    # HS
    df_HS_CyteGonia_stat = pd.DataFrame.from_records([
        ('HS', 'Cyte vs Gonia', 'Up', maxFDR, n_HS_UpCyteGonia_g),
        ('HS', 'Cyte vs Gonia', 'Not', maxFDR, n_HS_NotCyteGonia_g),
        ('HS', 'Cyte vs Gonia', 'Down', maxFDR, n_HS_DownCyteGonia_g)
    ], columns=columns)
    df_HS_CyteGonia_stat['%(G)'] = df_HS_CyteGonia_stat['Genes'] / df_HS_CyteGonia_stat['Genes'].sum()
    df_HS_CyteGonia_stat['Prot. Coding'] = [n_HS_UpCyteGonia_pcg, n_HS_NotCyteGonia_pcg, n_HS_DownCyteGonia_pcg]
    df_HS_CyteGonia_stat['%(PC)'] = df_HS_CyteGonia_stat['Prot. Coding'] / df_HS_CyteGonia_stat['Genes']

    # HS
    df_HS_CyteTid_stat = pd.DataFrame.from_records([
        ('HS', 'Cyte vs Tid', 'Up', maxFDR, n_HS_UpCyteTid_g),
        ('HS', 'Cyte vs Tid', 'Not', maxFDR, n_HS_NotCyteTid_g),
        ('HS', 'Cyte vs Tid', 'Down', maxFDR, n_HS_DownCyteTid_g)
    ], columns=columns)
    df_HS_CyteTid_stat['%(G)'] = df_HS_CyteTid_stat['Genes'] / df_HS_CyteTid_stat['Genes'].sum()
    df_HS_CyteTid_stat['Prot. Coding'] = [n_HS_UpCyteTid_pcg, n_HS_NotCyteTid_pcg, n_HS_DownCyteTid_pcg]
    df_HS_CyteTid_stat['%(PC)'] = df_HS_CyteTid_stat['Prot. Coding'] / df_HS_CyteTid_stat['Genes']

    # MM
    df_MM_CyteGonia_stat = pd.DataFrame.from_records([
        ('MM', 'Cyte vs Gonia', 'Up', maxFDR, n_MM_UpCyteGonia_g),
        ('MM', 'Cyte vs Gonia', 'Not', maxFDR, n_MM_NotCyteGonia_g),
        ('MM', 'Cyte vs Gonia', 'Down', maxFDR, n_MM_DownCyteGonia_g)
    ], columns=columns)
    df_MM_CyteGonia_stat['%(G)'] = df_MM_CyteGonia_stat['Genes'] / df_MM_CyteGonia_stat['Genes'].sum()
    df_MM_CyteGonia_stat['Prot. Coding'] = [n_MM_UpCyteGonia_pcg, n_MM_NotCyteGonia_pcg, n_MM_DownCyteGonia_pcg]
    df_MM_CyteGonia_stat['%(PC)'] = df_MM_CyteGonia_stat['Prot. Coding'] / df_MM_CyteGonia_stat['Genes']

    # MM
    df_MM_CyteTid_stat = pd.DataFrame.from_records([
        ('MM', 'Cyte vs Tid', 'Up', maxFDR, n_MM_UpCyteTid_g),
        ('MM', 'Cyte vs Tid', 'Not', maxFDR, n_MM_NotCyteTid_g),
        ('MM', 'Cyte vs Tid', 'Down', maxFDR, n_MM_DownCyteTid_g)
    ], columns=columns)
    df_MM_CyteTid_stat['%(G)'] = df_MM_CyteTid_stat['Genes'] / df_MM_CyteTid_stat['Genes'].sum()
    df_MM_CyteTid_stat['Prot. Coding'] = [n_MM_UpCyteTid_pcg, n_MM_NotCyteTid_pcg, n_MM_DownCyteTid_pcg]
    df_MM_CyteTid_stat['%(PC)'] = df_MM_CyteTid_stat['Prot. Coding'] / df_MM_CyteTid_stat['Genes']

    # DM
    df_DM_MiddleApical_stat = pd.DataFrame.from_records([
        ('DM', 'Middle vs Apical', 'Up', maxFDR, n_DM_UpMiddleApical_g),
        ('DM', 'Middle vs Apical', 'Not', maxFDR, n_DM_NotMiddleApical_g),
        ('DM', 'Middle vs Apical', 'Down', maxFDR, n_DM_DownMiddleApical_g)
    ], columns=columns)
    df_DM_MiddleApical_stat['%(G)'] = df_DM_MiddleApical_stat['Genes'] / df_DM_MiddleApical_stat['Genes'].sum()
    df_DM_MiddleApical_stat['Prot. Coding'] = [n_DM_UpMiddleApical_pcg, n_DM_NotMiddleApical_pcg, n_DM_DownMiddleApical_pcg]
    df_DM_MiddleApical_stat['%(PC)'] = df_DM_MiddleApical_stat['Prot. Coding'] / df_DM_MiddleApical_stat['Genes']

    # DM
    df_DM_MiddleBasal_stat = pd.DataFrame.from_records([
        ('DM', 'Middle vs Basal', 'Up', maxFDR, n_DM_UpMiddleBasal_g),
        ('DM', 'Middle vs Basal', 'Not', maxFDR, n_DM_NotMiddleBasal_g),
        ('DM', 'Middle vs Basal', 'Down', maxFDR, n_DM_DownMiddleBasal_g)
    ], columns=columns)
    df_DM_MiddleBasal_stat['%(G)'] = df_DM_MiddleBasal_stat['Genes'] / df_DM_MiddleBasal_stat['Genes'].sum()
    df_DM_MiddleBasal_stat['Prot. Coding'] = [n_DM_UpMiddleBasal_pcg, n_DM_NotMiddleBasal_pcg, n_DM_DownMiddleBasal_pcg]
    df_DM_MiddleBasal_stat['%(PC)'] = df_DM_MiddleBasal_stat['Prot. Coding'] / df_DM_MiddleBasal_stat['Genes']

    floatfmt = ['', '', '', '', '.2f', '', '.4f', '', '.4f']
    print(df2md(df_HS_CyteGonia_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_HS_CyteTid_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_MM_CyteGonia_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_MM_CyteTid_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_DM_MiddleApical_stat, floatfmt=floatfmt))
    print('')
    print(df2md(df_DM_MiddleBasal_stat, floatfmt=floatfmt))
    print('')