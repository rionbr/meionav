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
    # [D]rosophila [M]elanogaster (7227) - [A]liases
    #
    print('Mapping DM')
    # Query bioMart for Gene Name/Description
    ds_DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    df_DM_G = ds_DM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID')
    #
    rCSVFileUp = "../1-diff-gene-exp/results/DGE/DM/DM-DGE_Middle_vs_Apical.csv"
    rCSVFileDown = "../1-diff-gene-exp/results/DGE/DM/DM-DGE_Middle_vs_Basal.csv"
    df_DM_up = pd.read_csv(rCSVFileUp, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_DM_up.index.name = 'id_gene'
    df_DM_up.columns = [x + '_up' for x in df_DM_up.columns]
    df_DM_down = pd.read_csv(rCSVFileDown, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_DM_down.columns = [x + '_down' for x in df_DM_down.columns]
    df_DM_down.index.name = 'id_gene'

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_DM_up.index.to_list() + df_DM_down.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Up
    df_DM_up['id_string'] = df_SAg['id_string']
    df_DM_up['UpMiddle_vs_Apical'] = True
    # Down
    df_DM_down['id_string'] = df_SAg['id_string']
    df_DM_down['DownMiddle_vs_Basal'] = True

    # Merge Up/Down
    df_DM = pd.merge(df_DM_up, df_DM_down, how='outer', left_index=True, right_index=True)
    df_DM['id_string'] = df_DM.apply(combine_id_string_x_with_id_string_y, axis='columns')
    df_DM['gene'] = df_DM_G['Gene name']
    df_DM['biotype'] = df_DM_G['Gene type']
    df_DM[['UpMiddle_vs_Apical', 'DownMiddle_vs_Basal']] = df_DM[['UpMiddle_vs_Apical', 'DownMiddle_vs_Basal']].fillna(False)
    # Index Rows/Cols
    maskrows = (
        (df_DM['FDR_up'] <= 0.05) & (df_DM['logFC_up'].abs() >= minLogFC) & (df_DM['logFC_up'] >= 0) |
        (df_DM['FDR_down'] <= 0.05) & (df_DM['logFC_down'].abs() >= minLogFC) & (df_DM['logFC_down'] <= 0)
    )
    maskcols = [
        'id_string', 'gene', 'UpMiddle_vs_Apical', 'DownMiddle_vs_Basal',
        'logCPM_up', 'logFC_up', 'FDR_up',
        'logCPM_down', 'logFC_down', 'FDR_down',
        'biotype'
    ]
    df_DM = df_DM.loc[:, maskcols] # For Drosophila, we actually need all genes because of pipeline 'pooling'
    # To CSV
    df_DM.to_csv('results/DM-DE_genes.csv'.format(pipeline=pipeline))
    """
    #
    # [H]omo [S]apiens (9606) - [A]liases
    #
    """
    print('Mapping HS')
    # Query bioMart for Gene Name/Description
    ds_HS = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    df_HS_G = ds_HS.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')

    rCSVFileUp = "../1-diff-gene-exp/results/DGE/HS/HS-DGE_Cyte_vs_Gonia.csv"
    rCSVFileDown = "../1-diff-gene-exp/results/DGE/HS/HS-DGE_Cyte_vs_Tid.csv"
    df_HS_up = pd.read_csv(rCSVFileUp, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_HS_up.index.name = 'id_gene'
    df_HS_up.index = df_HS_up.index.map(lambda x: x.split('.')[0])
    df_HS_up.columns = [x + '_up' for x in df_HS_up.columns]
    df_HS_down = pd.read_csv(rCSVFileDown, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_HS_down.columns = [x + '_down' for x in df_HS_down.columns]
    df_HS_down.index.name = 'id_gene'
    df_HS_down.index = df_HS_down.index.map(lambda x: x.split('.')[0])
    
    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/9606/9606.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_HS_up.index.to_list() + df_HS_down.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Up
    df_HS_up['id_string'] = df_SAg['id_string']
    df_HS_up['UpCyte_vs_Gonia'] = True
    # Down
    df_HS_down['id_string'] = df_SAg['id_string']
    df_HS_down['DownCyte_vs_Tid'] = True

    # Merge Up/Down
    df_HS = pd.merge(df_HS_up, df_HS_down, how='outer', left_index=True, right_index=True)
    df_HS['id_string'] = df_HS.apply(combine_id_string_x_with_id_string_y, axis='columns')
    df_HS['gene'] = df_HS_G['Gene name']
    df_HS['biotype'] = df_HS_G['Gene type']
    df_HS[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']] = df_HS[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']].fillna(False)
    # Index Rows/Cols
    maskrows = (
        (df_HS['FDR_up'] <= 0.05) & (df_HS['logFC_up'].abs() >= minLogFC) & (df_HS['logFC_up'] >= 0) |
        (df_HS['FDR_down'] <= 0.05) & (df_HS['logFC_down'].abs() >= minLogFC) & (df_HS['logFC_down'] <= 0)
    )
    maskcols = [
        'id_string', 'gene', 'UpCyte_vs_Gonia', 'DownCyte_vs_Tid',
        'logCPM_up', 'logFC_up', 'FDR_up',
        'logCPM_down', 'logFC_down', 'FDR_down',
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

    rCSVFileUp = "../1-diff-gene-exp/results/DGE/MM/MM-DGE_Cyte_vs_Gonia.csv"
    rCSVFileDown = "../1-diff-gene-exp/results/DGE/MM/MM-DGE_Cyte_vs_Tid.csv"
    df_MM_up = pd.read_csv(rCSVFileUp, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_MM_up.index.name = 'id_gene'
    df_MM_up.index = df_MM_up.index.map(lambda x: x.split('.')[0])
    df_MM_up.columns = [x + '_up' for x in df_MM_up.columns]
    df_MM_down = pd.read_csv(rCSVFileDown, index_col=0).loc[:, ['logFC', 'logCPM', 'FDR']]
    df_MM_down.columns = [x + '_down' for x in df_MM_down.columns]
    df_MM_down.index.name = 'id_gene'
    df_MM_down.index = df_MM_down.index.map(lambda x: x.split('.')[0])

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/10090/10090.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_MM_up.index.to_list() + df_MM_down.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Up
    df_MM_up['id_string'] = df_SAg['id_string']
    df_MM_up['UpCyte_vs_Gonia'] = True
    # Down
    df_MM_down['id_string'] = df_SAg['id_string']
    df_MM_down['DownCyte_vs_Tid'] = True

    # Merge Up/Down
    df_MM = pd.merge(df_MM_up, df_MM_down, how='outer', left_index=True, right_index=True)
    df_MM['id_string'] = df_MM.apply(combine_id_string_x_with_id_string_y, axis='columns')
    df_MM['gene'] = df_MM_G['Gene name']
    df_MM['biotype'] = df_MM_G['Gene type']
    df_MM[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']] = df_MM[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']].fillna(False)
    # Index Rows/Cols
    maskrows = (
        (df_MM['FDR_up'] <= 0.05) & (df_MM['logFC_up'].abs() >= minLogFC) & (df_MM['logFC_up'] >= 0) |
        (df_MM['FDR_down'] <= 0.05) & (df_MM['logFC_down'].abs() >= minLogFC) & (df_MM['logFC_down'] <= 0)
    )
    maskcols = [
        'id_string', 'gene', 'UpCyte_vs_Gonia', 'DownCyte_vs_Tid',
        'logCPM_up', 'logFC_up', 'FDR_up',
        'logCPM_down', 'logFC_down', 'FDR_down',
        'biotype'
    ]
    df_MM = df_MM.loc[:, maskcols]
    # To CSV
    df_MM.to_csv('results/MM-DE_genes.csv'.format(pipeline=pipeline))

    print("Done.")
