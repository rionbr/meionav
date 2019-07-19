# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Generates the DM gene table that will be tested.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

if __name__ == '__main__':

    # Files to apply the map
    files = ['core_meiotic_genes_down.csv', 'core_meiotic_genes_up.csv'] # , 'core_meiotic_genes_updown.csv']

    # This is the 'df_SSG_3_up' object from [2]
    df_up = pd.read_csv('results/core_meiotic_genes_up.csv')
    df_down = pd.read_csv('results/core_meiotic_genes_down.csv')
    #df_updown = pd.read_csv('results/core_meiotic_genes_updown.csv')

    # Only what we are interested in
    cols = ['id_gene_HS', 'gene_HS', 'id_gene_DM', 'gene_DM', 'id_gene_MM', 'gene_MM']

    # For All three .csv files
    for ud, df in zip(['up', 'down'], [df_up, df_down]):
        df = df[cols]

        # Unpack DM genes
        unpckd = {'index': [], 'data': []}
        for (ids_gene_DM, genes_DM), row in df.set_index(['id_gene_DM', 'gene_DM']).iterrows():
            ids_gene_DM = ids_gene_DM.split(',')
            genes_DM = genes_DM.split(',')

            for (id_gene_DM, gene_DM) in zip(ids_gene_DM, genes_DM):
                unpckd['index'].append((id_gene_DM, gene_DM)) 
                unpckd['data'].append(row)

        dfU = pd.DataFrame(unpckd['data'], index=pd.MultiIndex.from_tuples(unpckd['index'], names=['id_gene_DM','gene_DM']))

        wCSVfile = 'results/core_DM_meiotic_genes_{:s}.csv'.format(ud)
        dfU.to_csv(wCSVfile)

    print('Done.')