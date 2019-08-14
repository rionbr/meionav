# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Processes original String-DB .zip.gz files.
#    Keeps only those ids we want.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import open_undefined_last_column_files, ensurePathExists

if __name__ == '__main__':

    # Compute results for both FDR <= 0.05 and FDR <= 0.01
    for maxFDR in [0.01, 0.05]:

        print("Computing for FDR <= {:f}".format(maxFDR))
        maxFDR_str = "FDR_{:s}".format(str(maxFDR).replace(".", "p"))

        # Load Species Files
        df_DM = pd.read_csv('results/DM/genes_DM.csv', index_col='id_string')
        rCSVFileHS = "results/HS/genes_HS-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        df_HS = pd.read_csv(rCSVFileHS, index_col='id_string')
        rCSVFileMM = "results/MM/genes_MM-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        df_MM = pd.read_csv(rCSVFileMM, index_col='id_string')

        # MapDict
        dfMap = pd.concat([df_HS, df_MM, df_DM], axis='index', sort=False)

        # Only Up Regulated
        string_HS = df_HS.index.tolist()
        string_MM = df_MM.index.tolist()
        string_DM = df_DM.index.tolist()

        string_HS_up = df_HS.loc[(df_HS['up'] == True), :].index.tolist()
        string_MM_up = df_MM.loc[(df_MM['up'] == True), :].index.tolist()
        string_DM_up = df_DM.loc[(df_DM['up'] == True), :].index.tolist()

        # Only Down Regulated
        string_HS_down = df_HS.loc[(df_HS['down'] == True), :].index.tolist()
        string_MM_down = df_MM.loc[(df_MM['down'] == True), :].index.tolist()
        string_DM_down = df_DM.loc[(df_DM['down'] == True), :].index.tolist()

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

        def selected_sperm_genes(genes,
                string_HS, string_MM, string_DM,
                string_HS_up, string_MM_up, string_DM_up,
                string_HS_down, string_MM_down, string_DM_down):
            # Separate by species, keeping only the sperm genes we are interested in
            spermgenes_HS = [gene for gene in genes if any([spermgene in gene for spermgene in string_HS])]
            spermgenes_MM = [gene for gene in genes if any([spermgene in gene for spermgene in string_MM])]
            spermgenes_DM = [gene for gene in genes if any([spermgene in gene for spermgene in string_DM])]

            # Mark as being Up or Down Regulated
            HS_up = [True if gene in string_HS_up else False for gene in spermgenes_HS]
            MM_up = [True if gene in string_MM_up else False for gene in spermgenes_MM]
            DM_up = [True if gene in string_DM_up else False for gene in spermgenes_DM]
            HS_down = [True if gene in string_HS_down else False for gene in spermgenes_HS]
            MM_down = [True if gene in string_MM_down else False for gene in spermgenes_MM]
            DM_down = [True if gene in string_DM_down else False for gene in spermgenes_DM]
            return pd.Series({
                'id_string_HS': spermgenes_HS, 'id_string_MM': spermgenes_MM, 'id_string_DM': spermgenes_DM,
                'HS_up': HS_up, 'MM_up': MM_up, 'DM_up': DM_up,
                'HS_down': HS_down, 'MM_down': MM_down, 'DM_down': DM_down,
            })

        print("> Separating by Species")
        df_SE = df_Egg.apply(separate_by_species)

        print("> Selecting Sperm Genes")
        df_SG = df_SE.apply(selected_sperm_genes,
            args=(
                string_HS, string_MM, string_DM,
                string_HS_up, string_MM_up, string_DM_up,
                string_HS_down, string_MM_down, string_DM_down
                )
            )

        # Map Gene ID
        def map_id_list(id_list, mapdict):
            """ Maps a list of id_string to a list of id_genes """
            return [mapdict[id_item] for id_item in id_list]

        # Map information
        print("> Mapping Additional Info")
        for species in ['HS', 'MM', 'DM']:
            # id_gene_HS, id_gene_MM, id_gene_DM
            df_SG['id_gene_' + species] = df_SG['id_string_' + species].apply(map_id_list, args=(dfMap['id_gene'].to_dict(), ))
            # gene_HS, gene_MM, gene_DM
            df_SG['gene_' + species] = df_SG['id_string_' + species].apply(map_id_list, args=(dfMap['gene'].to_dict(), ))
            # biotype_HS, biotype_MM, biotype_DM
            df_SG['biotype_' + species] = df_SG['id_string_' + species].apply(map_id_list, args=(dfMap['biotype'].to_dict(), ))
            # n_genes_HS, n_genes_MM, n_genes_DM
            df_SG['n_genes_' + species] = df_SG['id_gene_' + species].apply(len)
        df_SG['annotation'] = df_A['annotation']

        print("> Exporting")
        # Export this too
        wCSVFileSG = 'results/meiotic_genes/all_meiotic_genes-{maxFDR:s}.csv'.format(maxFDR=maxFDR_str)
        ensurePathExists(wCSVFileSG)
        df_SG.to_csv(wCSVFileSG)

        # Keep only genes with homolog in all three species
        df_SSG_3 = df_SG.loc[(df_SG[['n_genes_HS', 'n_genes_MM', 'n_genes_DM']] > 0).all(axis=1), :].copy()

        # Convert List to String
        columns = [
            'id_string_HS', 'id_string_MM', 'id_string_DM',
            'HS_up', 'MM_up', 'DM_up',
            'HS_down', 'MM_down', 'DM_down',
            'biotype_HS', 'biotype_MM', 'biotype_DM',
            'id_gene_HS', 'id_gene_MM', 'id_gene_DM',
            'gene_HS', 'gene_MM', 'gene_DM']
        for column in columns:
            df_SSG_3[column] = df_SSG_3[column].apply(lambda x: ",".join([str(y) for y in x]))

        wCSVFileSSG3 = 'results/core_meiotic_genes/core_meiotic_genes-{maxFDR:s}.csv'.format(maxFDR=maxFDR_str)
        ensurePathExists(wCSVFileSSG3)
        df_SSG_3.to_csv(wCSVFileSSG3)

    print('done.')
