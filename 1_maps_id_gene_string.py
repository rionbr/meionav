# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Processes original String-DB .zip.gz files.
#    Keeps only those genes that we want.
#
# NOTE: For some reason, "dmelanogaster_gene_ensembl" did not retrieve all gene names.
#       The following gene name was added manually after this script ran.
#       {"id_gene_DM"="FBgn0035707", "id_string_DM": "7227.FBpp0076685", "gene": "Rexo5"}
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import open_undefined_last_column_files
from pybiomart import Dataset

if __name__ == '__main__':

    def combine_id_string_x_with_id_string_y(row):
        if pd.isna(row['id_string_x']):
            return row['id_string_y']
        else:
            return row['id_string_x']

    #
    # [H]omo [S]apiens (9606) - [A]liases
    #
    id_genes_HS_up = pd.read_csv("data/HS_MeioticGenes_UpRegulated.csv", usecols=[0], squeeze=True, nrows=None).rename('id_gene_HS_up')
    id_genes_HS_down = pd.read_csv("data/HS_MeioticGenes_DownRegulated.csv", usecols=[0], squeeze=True, nrows=None).rename('id_gene_HS_down')

    df_HS_A = open_undefined_last_column_files('StringDB/9606/9606.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Map: id_gene <-> id_string
    df_HS_up = df_HS_A.loc[df_HS_A['alias'].isin(id_genes_HS_up), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_HS_up['up'] = True

    df_HS_down = df_HS_A.loc[df_HS_A['alias'].isin(id_genes_HS_down), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_HS_down['down'] = True

    df_HS = pd.merge(df_HS_up, df_HS_down, how='outer', left_index=True, right_index=True)
    df_HS['id_string'] = df_HS.apply(lambda r: r['id_string_x'] if pd.isna(r['id_string_y']) else r['id_string_y'], axis='columns')
    df_HS = df_HS[['id_string', 'up', 'down']]
    df_HS[['up', 'down']] = df_HS[['up', 'down']].fillna(False)

    # Query bioMart for Gene Name/Description
    ds_HS = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    df_HS_G = ds_HS.query(attributes=['ensembl_gene_id', 'external_gene_name', 'description']).set_index('Gene stable ID')

    df_HS['gene'] = df_HS_G['Gene name']
    df_HS['description'] = df_HS_G['Gene description']

    df_HS.to_csv('results/mapping/genes_HS.csv')

    #
    # [M]us [M]usculus (10090) - [A]liases
    #
    id_genes_MM_up = pd.read_csv("data/MM_MeioticGenes_UpRegulated.csv", usecols=[0], squeeze=True, nrows=None).rename('id_gene_MM_up')
    id_genes_MM_down = pd.read_csv("data/MM_MeioticGenes_DownRegulated.csv", usecols=[0], squeeze=True, nrows=None).rename('id_gene_MM_down')

    df_MM_A = open_undefined_last_column_files('StringDB/10090/10090.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Map: id_gene <-> id_string
    df_MM_up = df_MM_A.loc[df_MM_A['alias'].isin(id_genes_MM_up), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_MM_up['up'] = True

    df_MM_down = df_MM_A.loc[df_MM_A['alias'].isin(id_genes_MM_down), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_MM_down['down'] = True

    df_MM = pd.merge(df_MM_up, df_MM_down, how='outer', left_index=True, right_index=True)
    df_MM['id_string'] = df_MM.apply(lambda r: r['id_string_x'] if pd.isna(r['id_string_y']) else r['id_string_y'], axis='columns')
    df_MM = df_MM[['id_string', 'up', 'down']]
    df_MM[['up', 'down']] = df_MM[['up', 'down']].fillna(False)

    # Query bioMart for Gene Name/Description
    ds_MM = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
    df_MM_G = ds_MM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'description']).set_index('Gene stable ID')

    df_MM['gene'] = df_MM_G['Gene name']
    df_MM['description'] = df_MM_G['Gene description']

    df_MM.to_csv('results/mapping/genes_MM.csv')

    #
    # [D]rosophila [M]elanogaster (7227) - [A]liases
    #
    id_genes_DM_up = pd.read_csv("data/DM_MeioticGenes_UpRegulated.csv", usecols=[0], squeeze=True, nrows=None).rename('DS_idgene_up')
    id_genes_DM_down = pd.read_csv("data/DM_MeioticGenes_DownRegulated.csv", usecols=[0], squeeze=True, nrows=None).rename('DS_idgene_down')

    df_DM_A = open_undefined_last_column_files('StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Map: id_gene <-> id_string
    df_DM_up = df_DM_A.loc[df_DM_A['alias'].isin(id_genes_DM_up), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_DM_up['up'] = True

    df_DM_down = df_DM_A.loc[df_DM_A['alias'].isin(id_genes_DM_down), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_DM_down['down'] = True

    df_DM = pd.merge(df_DM_up, df_DM_down, how='outer', left_index=True, right_index=True)
    df_DM['id_string'] = df_DM.apply(lambda r: r['id_string_x'] if pd.isna(r['id_string_y']) else r['id_string_y'], axis='columns')
    df_DM = df_DM[['id_string', 'up', 'down']]
    df_DM[['up', 'down']] = df_DM[['up', 'down']].fillna(False)

    # Query bioMart for Gene Name/Description
    ds_DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    df_DM_G = ds_DM.query(attributes=['ensembl_gene_id', 'external_gene_name']).set_index('Gene stable ID')

    df_DM['gene'] = df_DM_G['Gene name']

    df_DM.to_csv('results/mapping/genes_DM.csv')

    print("Done.")
