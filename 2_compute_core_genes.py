# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Processes original String-DB .zip.gz files.
#    Keeps only those ids we want.
#
#
import re
import gzip
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import open_undefined_last_column_files

if __name__ == '__main__':

    # Load Species Files
    df_HS = pd.read_csv('results/mapping/genes_HS.csv', index_col='id_string')
    df_MM = pd.read_csv('results/mapping/genes_MM.csv', index_col='id_string')
    df_DM = pd.read_csv('results/mapping/genes_DM.csv', index_col='id_string')

    # MapDict
    dfMap = pd.concat([df_HS, df_MM, df_DM], axis='index', sort=False)

    # Only Up Regulated
    string_HS_up = df_HS.loc[(df_HS['up'] == True), :].index.tolist()
    string_MM_up = df_MM.loc[(df_MM['up'] == True), :].index.tolist()
    string_DM_up = df_DM.loc[(df_DM['up'] == True), :].index.tolist()

    # Only Down Regulated
    string_HS_down = df_HS.loc[(df_HS['down'] == True), :].index.tolist()
    string_MM_down = df_MM.loc[(df_MM['down'] == True), :].index.tolist()
    string_DM_down = df_DM.loc[(df_DM['down'] == True), :].index.tolist()

    # Both Up & Down Regulates
    #string_HS_updown = df_HS.loc[((df_HS['up'] == True) & (df_HS['down'] == True)), :].index.tolist()
    #string_MM_updown = df_MM.loc[((df_MM['up'] == True) & (df_MM['down'] == True)), :].index.tolist()
    #string_DM_updown = df_DM.loc[((df_DM['up'] == True) & (df_DM['down'] == True)), :].index.tolist()

    # Load EggNOG Annotation File
    df_A = pd.read_csv(
        'eggnog/33208_annotations.tsv',
        sep='\t',
        names=['species', 'id_eggnog', 'letter', 'annotation']).\
        set_index('id_eggnog')

    #
    # Metazoa (33208) EggNOG - [M]embers
    #
    df_Egg = open_undefined_last_column_files(
        "eggnog/33208_members.tsv.gz",
        n_fixed_cols=5,
        names=['family', 'id_eggnog', '_1', '_2', 'aliases', 'species'],
        nrows=None
    )

    # Only keep columns we need
    df_Egg = df_Egg.set_index('id_eggnog')['aliases']

    # List of species we are interested
    wanted_species = ['7227.', '9606.', '10090.']

    def separate_by_species(text):
        # Only keep genes from species we are interested in (lower the search space)
        return [gene for gene in text.split(',') if any([specie in gene for specie in wanted_species])]

    def selected_sperm_genes(genes, string_HS, string_MM, string_DM):
        # Separate by species, keeping only the sperm genes we are interested in
        spermgenes_HS = [gene for gene in genes if any([spermgene in gene for spermgene in string_HS])]
        spermgenes_MM = [gene for gene in genes if any([spermgene in gene for spermgene in string_MM])]
        spermgenes_DM = [gene for gene in genes if any([spermgene in gene for spermgene in string_DM])]
        #
        return pd.Series({'id_string_HS': spermgenes_HS, 'id_string_MM': spermgenes_MM, 'id_string_DM': spermgenes_DM})

    df_SE = df_Egg.apply(separate_by_species)
    df_SG_up = df_SE.apply(selected_sperm_genes, args=(string_HS_up, string_MM_up, string_DM_up))
    df_SG_down = df_SE.apply(selected_sperm_genes, args=(string_HS_down, string_MM_down, string_DM_down))
    #df_SG_updown = df_SE.apply(selected_sperm_genes, args=(string_HS_updown, string_MM_updown, string_DM_updown))

    # Map Gene ID
    def map_id_list(id_list, mapdict):
        """ Maps a list of id_string to a list of id_genes """
        return [mapdict[id_item] for id_item in id_list]

    # Map information to DataFrames
    for df in [df_SG_up, df_SG_down]: # , df_SG_updown]:
        # id_gene_<species>
        for species in ['HS', 'MM', 'DM']:
            # id_gene_HS, id_gene_MM, id_gene_DM
            df['id_gene_' + species] = df['id_string_' + species].apply(map_id_list, args=(dfMap['id_gene'].to_dict(), ))
            # gene_HS, gene_MM, gene_DM
            df['gene_' + species] = df['id_string_' + species].apply(map_id_list, args=(dfMap['gene'].to_dict(), ))
            # n_genes_HS, n_genes_MM, n_genes_DM
            df['n_genes_' + species] = df['id_gene_' + species].apply(len)
        df['annotation'] = df_A['annotation']

    """
    # Keep only genes with at least a homolog between 2 species
    dfSSG_up = df_SG_up.loc[
        (df_SG_up[['n_genes_DM', 'n_genes_HS']] > 0).all(axis=1) |
        (df_SG_up[['n_genes_DM', 'n_genes_MM']] > 0).all(axis=1) |
        (df_SG_up[['n_genes_HS', 'n_genes_MM']] > 0).all(axis=1), :].copy()
    dfSSG_down = df_SG_down.loc[
        (df_SG_down[['n_genes_DM', 'n_genes_HS']] > 0).all(axis=1) |
        (df_SG_down[['n_genes_DM', 'n_genes_MM']] > 0).all(axis=1) |
        (df_SG_down[['n_genes_HS', 'n_genes_MM']] > 0).all(axis=1), :].copy()
    """

    # Keep only genes with homolog in all three species
    df_SSG_3_up = df_SG_up.loc[
        (df_SG_up[['n_genes_HS', 'n_genes_MM', 'n_genes_DM']] > 0).all(axis=1), :].copy()
    df_SSG_3_down = df_SG_down.loc[
        (df_SG_down[['n_genes_HS', 'n_genes_MM', 'n_genes_DM']] > 0).all(axis=1), :].copy()
    #. df_SSG_3_updown = df_SG_updown.loc[
        # (df_SG_updown[['n_genes_HS', 'n_genes_MM', 'n_genes_DM']] > 0).all(axis=1), :].copy()

    # Convert List to String
    columns = [
        'id_string_HS', 'id_string_MM', 'id_string_DM',
        'id_gene_HS', 'id_gene_MM', 'id_gene_DM',
        'gene_HS', 'gene_MM', 'gene_DM']
    for column in columns:
        df_SSG_3_up[column] = df_SSG_3_up[column].apply(lambda x: ",".join([str(y) for y in x]))
        df_SSG_3_down[column] = df_SSG_3_down[column].apply(lambda x: ",".join([str(y) for y in x]))
        #df_SSG_3_updown[column] = df_SSG_3_updown[column].apply(lambda x: ",".join([str(y) for y in x]))

    df_SSG_3_up.to_csv('results/core_meiotic_genes_up.csv')
    df_SSG_3_down.to_csv('results/core_meiotic_genes_down.csv')
    # df_SSG_3_updown.to_csv('results/core_meiotic_genes_updown.csv')

    print('done.')
