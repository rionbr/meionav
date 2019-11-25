# coding=utf-8
# Author: Rion B Correia
# Date: Sept 11, 2019
#
# Description: Indexes meta-genes to select core meiotic genes.
#   Pipeline: Only mammal (HS & MM) conserved genes that Up/Down Regulated.
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists

if __name__ == '__main__':

    pipeline = 'mammals-conserved-FDRp05'
    maxFDR = 0.05
    minLogFC = math.log2(2)

    # Load Species Gene Files
    df_HS = pd.read_csv('results/HS-DE_genes.csv', index_col='id_gene')
    df_MM = pd.read_csv('results/MM-DE_genes.csv', index_col='id_gene')

    # Load Meta-Gene File
    df = pd.read_csv('results/meta_meiotic_genes.csv', index_col='id_eggnog')

    maskrows_HS = (
        (df_HS['Cyte_vs_Gonia'] == True) & (df_HS['FDR_CyteGonia'] <= maxFDR) & (df_HS['logFC_CyteGonia'].abs() >= minLogFC) & (df_HS['logFC_CyteGonia'] >= 0) |
        (df_HS['Tid_vs_Cyte'] == True) & (df_HS['FDR_TidCyte'] <= maxFDR) & (df_HS['logFC_TidCyte'].abs() >= minLogFC) & (df_HS['logFC_TidCyte'] <= 0)
    )
    maskrows_MM = (
        (df_MM['Cyte_vs_Gonia'] == True) & (df_MM['FDR_CyteGonia'] <= maxFDR) & (df_MM['logFC_CyteGonia'].abs() >= minLogFC) & (df_MM['logFC_CyteGonia'] >= 0) |
        (df_MM['Tid_vs_Cyte'] == True) & (df_MM['FDR_TidCyte'] <= maxFDR) & (df_MM['logFC_TidCyte'].abs() >= minLogFC) & (df_MM['logFC_TidCyte'] <= 0)
    )

    df_HS = df_HS.loc[maskrows_HS, :]
    df_MM = df_MM.loc[maskrows_MM, :]

    # List of id_strings
    string_HS = frozenset(np.hstack(df_HS['id_string'].values).tolist())
    string_MM = frozenset(np.hstack(df_MM['id_string'].values).tolist())

    def select_by_gene_and_species(text, keeplist):
        # Only keep genes we are interested in
        spermgenes = [i for i in text.split(',') if i in keeplist] if not pd.isna(text) else []
        return spermgenes if len(spermgenes) >= 1 else None

    print("> Selecting Sperm Genes")
    df['id_string_HS'] = df['id_string_HS'].apply(select_by_gene_and_species, args=(string_HS,))
    df['id_string_MM'] = df['id_string_MM'].apply(select_by_gene_and_species, args=(string_MM,))

    df = df.dropna(subset=['id_string_HS', 'id_string_MM'], how='any')

    # Map Gene ID
    def map_id_list(id_list, mapdict):
        """ Maps a list of id_string to a list of id_genes """
        return [mapdict[id_item] for id_item in id_list]

    """
    # Map information
    print("> Mapping Additional Info")
    for species in ['HS', 'MM', 'DM']:
        # id_gene_HS, id_gene_MM, id_gene_DM
        df['id_gene_' + species] = df['id_string_' + species].apply(map_id_list, args=(dfMap['id_gene'].to_dict(), ))
        # gene_HS, gene_MM, gene_DM
        df['gene_' + species] = df['id_string_' + species].apply(map_id_list, args=(dfMap['gene'].to_dict(), ))
        # biotype_HS, biotype_MM, biotype_DM
        df['biotype_' + species] = df['id_string_' + species].apply(map_id_list, args=(dfMap['biotype'].to_dict(), ))
        # n_genes_HS, n_genes_MM, n_genes_DM
        df['n_genes_' + species] = df['id_gene_' + species].apply(len)
    """
    # List of core id_strings
    core_string_HS = frozenset(np.hstack(df['id_string_HS'].values).tolist())
    core_string_MM = frozenset(np.hstack(df['id_string_MM'].values).tolist())

    df_HS = df_HS.loc[df_HS['id_string'].isin(core_string_HS), :]
    df_MM = df_MM.loc[df_MM['id_string'].isin(core_string_MM), :]

    # Convert List to String
    columns = ['id_string_HS', 'id_string_MM']
    for column in columns:
        df[column] = df[column].apply(lambda x: ",".join([str(y) for y in x]))

    columns = [
        'HS_CyteGonia', 'MM_CyteGonia',
        'HS_TidCyte', 'MM_TidCyte',
        'biotype_HS', 'biotype_MM',
        'id_gene_HS', 'id_gene_MM',
        'gene_HS', 'gene_MM']

    # Export
    print("> Exporting")
    wCSVFile = 'results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline)
    ensurePathExists(wCSVFile)
    df.to_csv(wCSVFile)
    
    # HS
    wCSVFileHS = 'results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline)
    ensurePathExists(wCSVFileHS)
    df_HS.to_csv(wCSVFileHS)

    # MM
    wCSVFileMM = 'results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline)
    ensurePathExists(wCSVFileMM)
    df_MM.to_csv(wCSVFileMM)

    print('done.')