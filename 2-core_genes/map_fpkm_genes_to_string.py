# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Maps FPKM genes to String-DB
#
import tarfile
import math
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import open_undefined_last_column_files
from pybiomart import Dataset


if __name__ == '__main__':

    #
    # [H]omo [S]apiens (9606) - [A]liases
    #
    print('Mapping HS')
    # Query bioMart for Gene Name/Description
    ds_HS = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
    df_HS_G = ds_HS.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')
    # Load FPKM
    tar_HS = tarfile.open('../1-diff-gene-exp/data/FPKM/HS/HS_Spermatocytes.tar.gz', "r:gz")
    df_HS = pd.concat(
        [
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep1.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep2.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep3.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep4.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep5.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep6.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep7.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep8.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep9.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep10.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep11.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_HS.extractfile('HS_Spermatocytes_rep12.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
        ],
        axis='index', sort=True)
    df_HS = df_HS.groupby(df_HS.index).mean()
    df_HS.index.name = 'id_gene'
    df_HS.index = df_HS.index.map(lambda x: x.split('.')[0])

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/9606/9606.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_HS.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Map
    df_HS['id_string'] = df_SAg['id_string']
    df_HS['gene'] = df_HS_G['Gene name']
    df_HS['Spermatocytes'] = True
    df_HS['biotype'] = df_HS_G['Gene type']
    # Index Rows/Cols
    maskcols = ['id_string', 'gene', 'Spermatocytes', 'FPKM', 'TPM', 'biotype']
    df_HS = df_HS.loc[:, maskcols]
    # To CSV
    df_HS.to_csv('results/HS-FPKM_genes.csv.gz')

    #
    # [M]us [M]usculus (10090) - [A]liases
    #
    print('Mapping MM')
    # Query bioMart for Gene Name/Description
    ds_MM = Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
    df_MM_G = ds_MM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID')
    # Load FPKM
    tar_MM = tarfile.open('../1-diff-gene-exp/data/FPKM/MM/MM_Spermatocytes.tar.gz', "r:gz")
    df_MM = pd.concat(
        [
            pd.read_csv(tar_MM.extractfile('MM_Spermatocytes_rep1.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_MM.extractfile('MM_Spermatocytes_rep2.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_MM.extractfile('MM_Spermatocytes_rep3.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
        ],
        axis='index', sort=True)
    df_MM = df_MM.groupby(df_MM.index).mean()
    df_MM.index.name = 'id_gene'

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/10090/10090.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])
    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_MM.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Map
    df_MM['id_string'] = df_SAg['id_string']
    df_MM['gene'] = df_MM_G['Gene name']
    df_MM['Spermatocytes'] = True
    df_MM['biotype'] = df_MM_G['Gene type']
    # Index Rows/Cols
    maskcols = ['id_string', 'gene', 'Spermatocytes', 'FPKM', 'TPM', 'biotype']
    df_MM = df_MM.loc[:, maskcols]
    # To CSV
    df_MM.to_csv('results/MM-FPKM_genes.csv.gz')

    #
    # [D]rosophila [M]elanogaster (7227) - [A]liases
    #
    print('Mapping DM')
    # Query bioMart for Gene Name/Description
    ds_DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    df_DM_G = ds_DM.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID')
    # Load FPKM
    tar_DM = tarfile.open('../1-diff-gene-exp/data/FPKM/DM/DM_Middle.tar.gz', "r:gz")
    df_DM = pd.concat(
        [
            pd.read_csv(tar_DM.extractfile('DM_Middle_rep1.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM']),
            pd.read_csv(tar_DM.extractfile('DM_Middle_rep2.txt'), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM'])
        ],
        axis='index', sort=True)
    df_DM = df_DM.groupby(df_DM.index).mean()
    df_DM.index.name = 'id_gene'

    # Map: id_gene <-> id_string
    df_SA = open_undefined_last_column_files('../StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"])

    # Parse String Data - Note some genes have multiple id_string, others have no match
    df_SA = df_SA.loc[df_SA['alias'].isin(df_DM.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
    df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
    # Map
    df_DM['id_string'] = df_SAg['id_string']
    df_DM['gene'] = df_DM_G['Gene name']
    df_DM['Middle'] = True
    df_DM['biotype'] = df_DM_G['Gene type']
    # Index Rows/Cols
    maskcols = ['id_string', 'gene', 'Middle', 'FPKM', 'TPM', 'biotype']
    df_DM = df_DM.loc[:, maskcols]
    # To CSV
    df_DM.to_csv('results/DM-FPKM_genes.csv.gz')

    print("Done.")
