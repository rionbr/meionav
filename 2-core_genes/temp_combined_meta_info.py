# coding=utf-8
# Author: Rion B Correia
# Date: Sept 19, 2019
#
# Description: Temporary file that joins the META file with gene id/name
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

    pipeline = 'all3-pooling-DM-FDRp05'
    # Load Files
    df = pd.read_csv('results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col='id_eggnog')
    df_HS = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col='id_string', usecols=['id_gene', 'id_string', 'gene'])
    df_MM = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col='id_string', usecols=['id_gene', 'id_string', 'gene'])
    df_DM = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col='id_string', usecols=['id_gene', 'id_string', 'gene'])

    def map_multiple_ids(x, d):
        x = x.split(',')
        return ','.join([d[i] for i in x])


    df['id_gene_HS'] = df['id_string_HS'].apply(map_multiple_ids, args=(df_HS['id_gene'].to_dict(),))
    df['id_gene_MM'] = df['id_string_MM'].apply(map_multiple_ids, args=(df_MM['id_gene'].to_dict(),))
    df['id_gene_DM'] = df['id_string_DM'].apply(map_multiple_ids, args=(df_DM['id_gene'].to_dict(),))

    df['gene_HS'] = df['id_string_HS'].apply(map_multiple_ids, args=(df_HS['gene'].to_dict(),))
    df['gene_MM'] = df['id_string_MM'].apply(map_multiple_ids, args=(df_MM['gene'].to_dict(),))
    df['gene_DM'] = df['id_string_DM'].apply(map_multiple_ids, args=(df_DM['gene'].to_dict(),))
    
    print("> Exporting")
    wCSVFile = 'results/{pipeline:s}/meta_meiotic_genes_4Paulo.csv'.format(pipeline=pipeline)
    ensurePathExists(wCSVFile)
    df.to_csv(wCSVFile)

    print('done.')
