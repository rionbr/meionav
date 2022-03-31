# coding=utf-8
# Author: Rion B Correia
# Date: Nov 12, 2021
#
# Description: Maps DGE genes to String-DB. Keeps only those genes that we want.
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

    #
    # [H]omo [S]apiens (9606) - [A]liases
    #

    print('Mapping HS')
    rCSVFile = "../../01-diff-gene-exp/results/rnf113/HS-DGE-rnf113_vs_control-SELECTED.csv"
    df_HS = pd.read_csv(rCSVFile, index_col=0)
    df_HS.index.name = 'id_gene'

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../../data/StringDB/9606/9606.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_HS.index), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    #
    df_HS['id_string'] = df_SAg['id_string']

    # To CSV
    df_HS.to_csv('results/HS-DGE-rnf113_vs_control-SELECTED-string.csv')

    #
    # [M]us [M]usculus (10090) - [A]liases
    #
    """
    print('Mapping MM')
    rCSVFile = "../01-diff-gene-exp/results/MM/MM-DGE_Cyte_vs_Gonia.csv"
    df_MM = pd.read_csv(rCSVFile, index_col=0)
    df_MM.index.name = 'id_gene'

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../data/StringDB/10090/10090.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_MM.index), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    #
    df_MM['id_string'] = df_SAg['id_string']

    # To CSV
    df_MM.to_csv('results/MM-DGE-rnf113_vs_control-SELECTED-string.csv')
    """
    #
    # [D]rosophila [M]elanogaster (7227) - [A]liases
    #
    print('Mapping DM')
    #
    rCSVFile = "../../01-diff-gene-exp/results/mdlc/DM-DGE-mdlc_vs_control-SELECTED.csv"
    df_DM = pd.read_csv(rCSVFile, index_col=0)
    df_DM.index.name = 'id_gene'

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../../data/StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])

    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_DM.index), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    #
    df_DM['id_string'] = df_SAg['id_string']

    # To CSV
    df_DM.to_csv('results/DM-DGE-mdlc_vs_control-SELECTED-string.csv')

    print("Done.")
