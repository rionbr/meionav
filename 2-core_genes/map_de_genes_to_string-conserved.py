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

    pipeline = 'conserved'

    # Compute results for both FDR <= 0.05 and FDR <= 0.01
    for maxFDR in [0.01, 0.05]:

        print("- FDR <= {:f}".format(maxFDR))
        maxFDR_str = "FDR_{:s}".format(str(maxFDR).replace(".", "p"))
        #
        # [D]rosophila [M]elanogaster (7227) - [A]liases
        #
        print('Mapping DM')
        rCSVFileUp = "../1-diff-gene-exp/results/DE-{pipeline:s}/DM/DM_DE_UpMiddle_vs_Apical-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        rCSVFileDown = "../1-diff-gene-exp/results/DE-{pipeline:s}/DM/DM_DE_DownMiddle_vs_Basal-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        id_genes_DM_up = pd.read_csv(rCSVFileUp, usecols=[0], squeeze=True, nrows=None).rename('id_gene_DM_up')
        id_genes_DM_down = pd.read_csv(rCSVFileDown, usecols=[0], squeeze=True, nrows=None).rename('id_gene_DM_down')

        df_DM_A = open_undefined_last_column_files('../StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
        # Map: id_gene <-> id_string
        df_DM_up = df_DM_A.loc[df_DM_A['alias'].isin(id_genes_DM_up), ["alias", "id_string"]].\
            rename(columns={"alias": "id_gene"}).\
            set_index('id_gene', drop=True)
        df_DM_up['UpMiddle_vs_Apical'] = True

        df_DM_down = df_DM_A.loc[df_DM_A['alias'].isin(id_genes_DM_down), ["alias", "id_string"]].\
            rename(columns={"alias": "id_gene"}).\
            set_index('id_gene', drop=True)
        df_DM_down['DownMiddle_vs_Basal'] = True

        df_DM = pd.merge(df_DM_up, df_DM_down, how='outer', left_index=True, right_index=True)
        df_DM['id_string'] = df_DM.apply(lambda r: r['id_string_x'] if pd.isna(r['id_string_y']) else r['id_string_y'], axis='columns')
        df_DM = df_DM[['id_string', 'UpMiddle_vs_Apical', 'DownMiddle_vs_Basal']]
        df_DM[['UpMiddle_vs_Apical', 'DownMiddle_vs_Basal']] = df_DM[['UpMiddle_vs_Apical', 'DownMiddle_vs_Basal']].fillna(False)

        # Query bioMart for Gene Name/Description
        ds_DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
        df_DM_G = ds_DM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID')

        # Manual Inclusion
        # df_DM_G.loc['FBgn0035707', 'Gene name'] = 'Rexo5'
        # df_DM_G.loc['FBgn0038035', 'Gene name'] = 'CG17227 (DNAlig3)'
        # df_DM_G.loc['FBgn0038341', 'Gene name'] = 'CG14869 (AdamTS-A)'

        df_DM['gene'] = df_DM_G['Gene name']
        df_DM['biotype'] = df_DM_G['Gene type']

        wCSVFileDM = 'results/{pipeline:s}/DM/genes_DM-{maxFDR:s}.csv'.format(pipeline=pipeline, maxFDR=maxFDR_str)
        ensurePathExists(wCSVFileDM)
        df_DM.to_csv(wCSVFileDM)

        #
        # [H]omo [S]apiens (9606) - [A]liases
        #
        print('Mapping HS')
        rCSVFileUp = "../1-diff-gene-exp/results/DE-{pipeline:s}/HS/HS_DE_UpCyte_vs_Gonia-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        rCSVFileDown = "../1-diff-gene-exp/results/DE-{pipeline:s}/HS/HS_DE_DownCyte_vs_Tid-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        id_genes_HS_up = pd.read_csv(rCSVFileUp, usecols=[0], squeeze=True, nrows=None).rename('id_gene_HS_up')
        id_genes_HS_down = pd.read_csv(rCSVFileDown, usecols=[0], squeeze=True, nrows=None).rename('id_gene_HS_down')

        # Remove gene version
        id_genes_HS_up = id_genes_HS_up.map(lambda x: x.split('.')[0])
        id_genes_HS_down = id_genes_HS_down.map(lambda x: x.split('.')[0])

        df_HS_A = open_undefined_last_column_files('../StringDB/9606/9606.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
        # Map: id_gene <-> id_string
        df_HS_up = df_HS_A.loc[df_HS_A['alias'].isin(id_genes_HS_up), ["alias", "id_string"]].\
            rename(columns={"alias": "id_gene"}).\
            set_index('id_gene', drop=True)
        df_HS_up['UpCyte_vs_Gonia'] = True

        df_HS_down = df_HS_A.loc[df_HS_A['alias'].isin(id_genes_HS_down), ["alias", "id_string"]].\
            rename(columns={"alias": "id_gene"}).\
            set_index('id_gene', drop=True)
        df_HS_down['DownCyte_vs_Tid'] = True

        df_HS = pd.merge(df_HS_up, df_HS_down, how='outer', left_index=True, right_index=True)

        df_HS['id_string'] = df_HS.apply(lambda r: r['id_string_x'] if pd.isna(r['id_string_y']) else r['id_string_y'], axis='columns')
        df_HS = df_HS[['id_string', 'UpCyte_vs_Gonia', 'DownCyte_vs_Tid']]
        df_HS[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']] = df_HS[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']].fillna(False)

        # Query bioMart for Gene Name/Description
        ds_HS = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
        df_HS_G = ds_HS.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')

        df_HS['gene'] = df_HS_G['Gene name']
        df_HS['biotype'] = df_HS_G['Gene type']
        df_HS['description'] = df_HS_G['Gene description']

        wCSVFileHS = "results/{pipeline:s}/HS/genes_HS-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        ensurePathExists(wCSVFileHS)
        df_HS.to_csv(wCSVFileHS)

        #
        # [M]us [M]usculus (10090) - [A]liases
        #
        print('Mapping MM')
        rCSVFileUp = "../1-diff-gene-exp/results/DE-{pipeline:s}/MM/MM_DE_UpCyte_vs_Gonia-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        rCSVFileDown = "../1-diff-gene-exp/results/DE-{pipeline:s}/MM/MM_DE_DownCyte_vs_Tid-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        id_genes_MM_up = pd.read_csv(rCSVFileUp, usecols=[0], squeeze=True, nrows=None).rename('id_gene_MM_up')
        id_genes_MM_down = pd.read_csv(rCSVFileDown, usecols=[0], squeeze=True, nrows=None).rename('id_gene_MM_down')

        df_MM_A = open_undefined_last_column_files('../StringDB/10090/10090.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
        # Map: id_gene <-> id_string
        df_MM_up = df_MM_A.loc[df_MM_A['alias'].isin(id_genes_MM_up), ["alias", "id_string"]].\
            rename(columns={"alias": "id_gene"}).\
            set_index('id_gene', drop=True)
        df_MM_up['UpCyte_vs_Gonia'] = True

        df_MM_down = df_MM_A.loc[df_MM_A['alias'].isin(id_genes_MM_down), ["alias", "id_string"]].\
            rename(columns={"alias": "id_gene"}).\
            set_index('id_gene', drop=True)
        df_MM_down['DownCyte_vs_Tid'] = True

        df_MM = pd.merge(df_MM_up, df_MM_down, how='outer', left_index=True, right_index=True)
        df_MM['id_string'] = df_MM.apply(lambda r: r['id_string_x'] if pd.isna(r['id_string_y']) else r['id_string_y'], axis='columns')
        df_MM = df_MM[['id_string', 'UpCyte_vs_Gonia', 'DownCyte_vs_Tid']]
        df_MM[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']] = df_MM[['UpCyte_vs_Gonia', 'DownCyte_vs_Tid']].fillna(False)

        # Query bioMart for Gene Name/Description
        ds_MM = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
        df_MM_G = ds_MM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')

        df_MM['gene'] = df_MM_G['Gene name']
        df_MM['biotype'] = df_MM_G['Gene type']
        df_MM['description'] = df_MM_G['Gene description']

        wCSVFileMM = "results/{pipeline:s}/MM/genes_MM-{maxFDR:s}.csv".format(pipeline=pipeline, maxFDR=maxFDR_str)
        ensurePathExists(wCSVFileMM)
        df_MM.to_csv(wCSVFileMM)

    print("Done.")
