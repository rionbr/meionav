# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Generates the DM gene table that will be tested.
#   Also generates the possible other species genes to tets.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists


if __name__ == '__main__':

    pipeline = 'pooling'
    ldfDM = []
    ldfHS = []
    ldfMM = []

    for maxFDR in [0.05, 0.01]:

        print("- FDR <= {:f}".format(maxFDR))
        maxFDR_str = "FDR_{:s}".format(str(maxFDR).replace(".", "p"))

        rCSVFile = 'results/{pipeline:s}/core_meiotic_genes/core_meiotic_genes-{maxFDR:s}.csv'.format(pipeline=pipeline, maxFDR=maxFDR_str)
        df = pd.read_csv(rCSVFile)

        # Only what we are interested in
        cols = [
            'id_gene_HS', 'gene_HS', 'biotype_HS', 'HS_up', 'HS_down',
            'id_gene_DM', 'gene_DM', 'biotype_DM', 'DM_ma', 'DM_mb',  # changes here
            'id_gene_MM', 'gene_MM', 'biotype_MM', 'MM_up', 'MM_down'
        ]

        # ReArrange Columns
        df = df[cols]

        #
        # Unpack DM genes
        #
        unpckd = {'index': [], 'data': []}
        for (ids_gene_DM, genes_DM, biotypes_DM, DM_mas, DM_mbs), row in df.set_index(['id_gene_DM', 'gene_DM', 'biotype_DM', 'DM_ma', 'DM_mb']).iterrows():
            ids_gene_DM = ids_gene_DM.split(',')
            genes_DM = genes_DM.split(',')
            biotypes_DM = biotypes_DM.split(',')
            DM_mas = DM_mas.split(',')
            DM_mbs = DM_mbs.split(',')

            for (id_gene_DM, gene_DM, biotype_DM, DM_ma, DM_mb) in zip(ids_gene_DM, genes_DM, biotypes_DM, DM_mas, DM_mbs):
                unpckd['index'].append((id_gene_DM, gene_DM, biotype_DM, DM_ma, DM_mb))
                unpckd['data'].append(row)

        dfU = pd.DataFrame(unpckd['data'], index=pd.MultiIndex.from_tuples(unpckd['index'], names=['id_gene_DM', 'gene_DM', 'biotype_DM', 'DM_ma', 'DM_mb']))

        # Sometimes two different proteins matched in EggNOG map to the same gene. This just removes the duplicates
        dfU = dfU.reset_index(level=['gene_DM', 'DM_ma', 'DM_mb', 'biotype_DM']).drop_duplicates(keep='first')
        ldfDM.append(dfU)

        #
        # Unpack HS genes
        #
        unpckd = {'index': [], 'data': []}
        for (ids_gene_HS, genes_HS, biotypes_HS, HS_ups, HS_downs), row in df.set_index(['id_gene_HS', 'gene_HS', 'biotype_HS', 'HS_up', 'HS_down']).iterrows():
            ids_gene_HS = ids_gene_HS.split(',')
            genes_HS = genes_HS.split(',')
            biotypes_HS = biotypes_HS.split(',')
            HS_ups = HS_ups.split(',')
            HS_downs = HS_downs.split(',')

            for (id_gene_HS, gene_HS, biotype_HS, HS_up, HS_down) in zip(ids_gene_HS, genes_HS, biotypes_HS, HS_ups, HS_downs):
                unpckd['index'].append((id_gene_HS, gene_HS, biotype_HS, HS_up, HS_down))
                unpckd['data'].append(row)

        dfU = pd.DataFrame(unpckd['data'], index=pd.MultiIndex.from_tuples(unpckd['index'], names=['id_gene_HS', 'gene_HS', 'biotype_HS', 'HS_up', 'HS_down']))

        # Sometimes two different proteins matched in EggNOG map to the same gene. This just removes the duplicates
        dfU = dfU.reset_index(level=['gene_HS', 'HS_up', 'HS_down', 'biotype_HS']).drop_duplicates(keep='first')
        ldfHS.append(dfU)

        #
        # Unpack MM genes
        #
        unpckd = {'index': [], 'data': []}
        for (ids_gene_MM, genes_MM, biotypes_MM, MM_ups, MM_downs), row in df.set_index(['id_gene_MM', 'gene_MM', 'biotype_MM', 'MM_up', 'MM_down']).iterrows():
            ids_gene_MM = ids_gene_MM.split(',')
            genes_MM = genes_MM.split(',')
            biotypes_MM = biotypes_MM.split(',')
            MM_ups = MM_ups.split(',')
            MM_downs = MM_downs.split(',')

            for (id_gene_MM, gene_MM, biotype_MM, MM_up, MM_down) in zip(ids_gene_MM, genes_MM, biotypes_MM, MM_ups, MM_downs):
                unpckd['index'].append((id_gene_MM, gene_MM, biotype_MM, MM_up, MM_down))
                unpckd['data'].append(row)

        dfU = pd.DataFrame(unpckd['data'], index=pd.MultiIndex.from_tuples(unpckd['index'], names=['id_gene_MM', 'gene_MM', 'biotype_MM', 'MM_up', 'MM_down']))

        # Sometimes two different proteins matched in EggNOG map to the same gene. This just removes the duplicates
        dfU = dfU.reset_index(level=['gene_MM', 'MM_up', 'MM_down', 'biotype_MM']).drop_duplicates(keep='first')
        ldfMM.append(dfU)


    #
    # Generate one table that contains it all
    #
    dfDM = ldfDM[0]
    dfDM1 = ldfDM[1]
    dfDM['FDR<=0.05'] = True
    dfDM['FDR<=0.01'] = dfDM.index.isin(dfDM1.index)
    dfDM.to_csv('results/{pipeline:s}/DM/core_DM_meiotic_genes.csv'.format(pipeline=pipeline))

    dfHS = ldfHS[0]
    dfHS1 = ldfHS[1]
    dfHS['FDR<=0.05'] = True
    dfHS['FDR<=0.01'] = dfHS.index.isin(dfHS1.index)
    dfHS.to_csv('results/{pipeline:s}/HS/core_HS_meiotic_genes.csv'.format(pipeline=pipeline))

    dfMM = ldfMM[0]
    dfMM1 = ldfMM[1]
    dfMM['FDR<=0.05'] = True
    dfMM['FDR<=0.01'] = dfMM.index.isin(dfMM1.index)
    dfMM.to_csv('results/{pipeline:s}/MM/core_MM_meiotic_genes.csv'.format(pipeline=pipeline))
