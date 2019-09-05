# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Maps DE genes to String-DB. Keeps only those genes that we want.
#
# NOTE: For some reason, "dmelanogaster_gene_ensembl" did not retrieve all gene names. Some were manually added at the end.
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import open_undefined_last_column_files, ensurePathExists
from pybiomart import Dataset


def combine_id_string_x_with_id_string_y(row):
    if pd.isna(row['id_string_x']):
        return row['id_string_y']
    else:
        return row['id_string_x']


if __name__ == '__main__':

    pipeline = 'pooling'

    #
    # [D]rosophila [M]elanogaster (7227) - [A]liases
    #
    print('Mapping DM')
    rCSVFileMiddleApical = "../1-diff-gene-exp/results/DE-{pipeline:s}/DM/DM_DE_Middle_vs_Apical.csv".format(pipeline=pipeline)
    rCSVFileMiddleBasal = "../1-diff-gene-exp/results/DE-{pipeline:s}/DM/DM_DE_Middle_vs_Basal.csv".format(pipeline=pipeline)
    id_genes_DM_MA = pd.read_csv(rCSVFileMiddleApical, usecols=[0], squeeze=True, nrows=None).rename('id_gene_DM')
    id_genes_DM_MB = pd.read_csv(rCSVFileMiddleBasal, usecols=[0], squeeze=True, nrows=None).rename('id_gene_DM')

    id_genes_DM = pd.concat([id_genes_DM_MA, id_genes_DM_MB], axis=0).drop_duplicates().reset_index(drop=True)
    
    df_DM_A = open_undefined_last_column_files('../StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Map: id_gene <-> id_string
    df_DM_up = df_DM_A.loc[df_DM_A['alias'].isin(id_genes_DM_MA), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_DM_up['Middle_vs_Apical'] = True

    df_DM_down = df_DM_A.loc[df_DM_A['alias'].isin(id_genes_DM_MB), ["alias", "id_string"]].\
        rename(columns={"alias": "id_gene"}).\
        set_index('id_gene', drop=True)
    df_DM_down['Middle_vs_Basal'] = True

    df_DM = pd.merge(df_DM_up, df_DM_down, how='outer', left_index=True, right_index=True)
    df_DM['id_string'] = df_DM.apply(lambda r: r['id_string_x'] if pd.isna(r['id_string_y']) else r['id_string_y'], axis='columns')
    df_DM = df_DM[['id_string', 'Middle_vs_Apical', 'Middle_vs_Basal']]
    df_DM[['Middle_vs_Apical', 'Middle_vs_Basal']] = df_DM[['Middle_vs_Apical', 'Middle_vs_Basal']].fillna(False)

    # Query bioMart for Gene Name/Description
    ds_DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    df_DM_G = ds_DM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID')

    # Manual Inclusion
    # df_DM_G.loc['FBgn0035707', 'Gene name'] = 'Rexo5'
    # df_DM_G.loc['FBgn0038035', 'Gene name'] = 'CG17227 (DNAlig3)'
    # df_DM_G.loc['FBgn0038341', 'Gene name'] = 'CG14869 (AdamTS-A)'

    df_DM['gene'] = df_DM_G['Gene name']
    df_DM['biotype'] = df_DM_G['Gene type']

    wCSVFileDM = 'results/{pipeline:s}/DM/genes_DM.csv'.format(pipeline=pipeline)
    ensurePathExists(wCSVFileDM)
    df_DM.to_csv(wCSVFileDM)

    print("Done.")
