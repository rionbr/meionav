# coding=utf-8
# Author: Rion B Correia
# Date: Ago 14, 2019
#
# Description: Extracts the list of genes that were miss-screened based on a previous pipeline.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

if __name__ == '__main__':

    df = pd.read_csv('../2-core_genes/results/DM/core_DM_meiotic_genes.csv').set_index('id_gene_DM')
    
    # Screened Data (From Experimental Analysis)
    dfS = pd.read_csv('data/core_DM_screened.csv').set_index('id_gene_DM')

    dfMISS = dfS.loc[~dfS.index.isin(df.index),:].index.tolist()
    df = df.loc[notin,:]
