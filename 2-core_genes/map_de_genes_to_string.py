# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Maps DE genes to String-DB. Keeps only those genes that we want.
#
# NOTE: For some reason, "dmelanogaster_gene_ensembl" did not retrieve all gene names. Some were manually added at the end.
#
import math
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import open_undefined_last_column_files, ensurePathExists
from pybiomart import Dataset


def combine_id_string_x_with_id_string_y(r):
    x = r['id_string_x']
    y = r['id_string_y']
    if isinstance(x, list):
        return x
    elif not pd.isna(x):
        return x
    else:
        return y

if __name__ == '__main__':

    pipeline = 'conserved'
    minLogFC = math.log2(2)


    #
    # [H]omo [S]apiens (9606) - [A]liases
    #
    print('Mapping HS')
    # Query bioMart for Gene Name/Description
    ds_HS = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    df_HS_G = ds_HS.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')

    rCSVFileCG = "../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Gonia.csv"
    rCSVFileCT = "../1-diff-gene-exp/results/HS/HS-DGE_Cyte_vs_Tid.csv"
    df_HS_CG = pd.read_csv(rCSVFileCG, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_HS_CG.index.name = 'id_gene'
    df_HS_CG.index = df_HS_CG.index.map(lambda x: x.split('.')[0])
    df_HS_CG.columns = [x + '_CyteGonia' for x in df_HS_CG.columns]
    df_HS_CT = pd.read_csv(rCSVFileCT, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_HS_CT.columns = [x + '_CyteTid' for x in df_HS_CT.columns]
    df_HS_CT.index.name = 'id_gene'
    df_HS_CT.index = df_HS_CT.index.map(lambda x: x.split('.')[0])
    
    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/9606/9606.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_HS_CG.index.to_list() + df_HS_CT.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Up
    df_HS_CG['id_string'] = df_SAg['id_string']
    df_HS_CG['Cyte_vs_Gonia'] = True
    # Down
    df_HS_CT['id_string'] = df_SAg['id_string']
    df_HS_CT['Cyte_vs_Tid'] = True

    # Merge Up/Down
    df_HS = pd.merge(df_HS_CG, df_HS_CT, how='outer', left_index=True, right_index=True)
    df_HS['id_string'] = df_HS.apply(combine_id_string_x_with_id_string_y, axis='columns')
    df_HS['gene'] = df_HS_G['Gene name']
    df_HS['biotype'] = df_HS_G['Gene type']
    df_HS[['Cyte_vs_Gonia', 'Cyte_vs_Tid']] = df_HS[['Cyte_vs_Gonia', 'Cyte_vs_Tid']].fillna(False)
    # Index Rows/Cols
    maskcols = [
        'id_string', 'gene', 'Cyte_vs_Gonia', 'Cyte_vs_Tid',
        'logCPM_CyteGonia', 'logFC_CyteGonia', 'FDR_CyteGonia',
        'logCPM_CyteTid', 'logFC_CyteTid', 'FDR_CyteTid',
        'biotype'
    ]
    df_HS = df_HS.loc[:, maskcols]
    # To CSV
    df_HS.to_csv('results/HS-DE_genes.csv'.format(pipeline=pipeline))
    
    #
    # [M]us [M]usculus (10090) - [A]liases
    #
    print('Mapping MM')
    # Query bioMart for Gene Name/Description
    ds_MM = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
    df_MM_G = ds_MM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')

    rCSVFileCG = "../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv"
    rCSVFileCT = "../1-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Tid.csv"
    df_MM_CG = pd.read_csv(rCSVFileCG, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_MM_CG.index.name = 'id_gene'
    df_MM_CG.index = df_MM_CG.index.map(lambda x: x.split('.')[0])
    df_MM_CG.columns = [x + '_CyteGonia' for x in df_MM_CG.columns]
    df_MM_CT = pd.read_csv(rCSVFileCT, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_MM_CT.columns = [x + '_CyteTid' for x in df_MM_CT.columns]
    df_MM_CT.index.name = 'id_gene'
    df_MM_CT.index = df_MM_CT.index.map(lambda x: x.split('.')[0])

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/10090/10090.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_MM_CG.index.to_list() + df_MM_CT.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Up
    df_MM_CG['id_string'] = df_SAg['id_string']
    df_MM_CG['Cyte_vs_Gonia'] = True
    # Down
    df_MM_CT['id_string'] = df_SAg['id_string']
    df_MM_CT['Cyte_vs_Tid'] = True

    # Merge Up/Down
    df_MM = pd.merge(df_MM_CG, df_MM_CT, how='outer', left_index=True, right_index=True)
    df_MM['id_string'] = df_MM.apply(combine_id_string_x_with_id_string_y, axis='columns')
    df_MM['gene'] = df_MM_G['Gene name']
    df_MM['biotype'] = df_MM_G['Gene type']
    df_MM[['Cyte_vs_Gonia', 'Cyte_vs_Tid']] = df_MM[['Cyte_vs_Gonia', 'Cyte_vs_Tid']].fillna(False)
    # Index Rows/Cols
    maskcols = [
        'id_string', 'gene', 'Cyte_vs_Gonia', 'Cyte_vs_Tid',
        'logCPM_CyteGonia', 'logFC_CyteGonia', 'FDR_CyteGonia',
        'logCPM_CyteTid', 'logFC_CyteTid', 'FDR_CyteTid',
        'biotype'
    ]
    df_MM = df_MM.loc[:, maskcols]
    # To CSV
    df_MM.to_csv('results/MM-DE_genes.csv'.format(pipeline=pipeline))

    #
    # [D]rosophila [M]elanogaster (7227) - [A]liases
    #
    print('Mapping DM')
    # Query bioMart for Gene Name/Description
    ds_DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    df_DM_G = ds_DM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID')
    #
    rCSVFileMA = "../1-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Apical.csv"
    rCSVFileMB = "../1-diff-gene-exp/results/DM/DM-DGE_Middle_vs_Basal.csv"
    df_DM_MA = pd.read_csv(rCSVFileMA, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_DM_MA.index.name = 'id_gene'
    df_DM_MA.columns = [x + '_MiddleApical' for x in df_DM_MA.columns]
    df_DM_MB = pd.read_csv(rCSVFileMB, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_DM_MB.columns = [x + '_MiddleBasal' for x in df_DM_MB.columns]
    df_DM_MB.index.name = 'id_gene'

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_DM_MA.index.to_list() + df_DM_MB.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Up
    df_DM_MA['id_string'] = df_SAg['id_string']
    df_DM_MA['Middle_vs_Apical'] = True
    # Down
    df_DM_MB['id_string'] = df_SAg['id_string']
    df_DM_MB['Middle_vs_Basal'] = True

    # Merge Up/Down
    df_DM = pd.merge(df_DM_MA, df_DM_MB, how='outer', left_index=True, right_index=True)
    df_DM['id_string'] = df_DM.apply(combine_id_string_x_with_id_string_y, axis='columns')
    df_DM['gene'] = df_DM_G['Gene name']
    df_DM['biotype'] = df_DM_G['Gene type']
    df_DM[['Middle_vs_Apical', 'Middle_vs_Basal']] = df_DM[['Middle_vs_Apical', 'Middle_vs_Basal']].fillna(False)
    # Index Rows/Cols
    maskcols = [
        'id_string', 'gene', 'Middle_vs_Apical', 'Middle_vs_Basal',
        'logCPM_MiddleApical', 'logFC_MiddleApical', 'FDR_MiddleApical',
        'logCPM_MiddleBasal', 'logFC_MiddleBasal', 'FDR_MiddleBasal',
        'biotype'
    ]
    df_DM = df_DM.loc[:, maskcols] # For Drosophila, we actually need all genes because of pipeline 'pooling'
    # To CSV
    df_DM.to_csv('results/DM-DE_genes.csv'.format(pipeline=pipeline))
    
    print("Done.")
