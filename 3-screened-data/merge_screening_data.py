# coding=utf-8
# Author: Rion B Correia
# Date: Jul 17, 2019
#
# Description: Merges Selected Gens with Screening Data
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

if __name__ == '__main__':

    # This is the 'dfU' object from [3]
    df = pd.read_csv('../2-core_genes/results/DM/core_DM_meiotic_genes.csv').set_index('id_gene_DM')
    
    # Screened Data (From Experimental Analysis)
    dfS = pd.read_csv('data/core_DM_screened_2019-08.csv').set_index('id_gene_DM')
    print(dfS.head())
    cols1 = [
        'StockNumber', 'FlyBase Genotype'
    ]
    cols2 = [
        'FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched',
    ]
    cols3 = [
        'Phenotype or literature', 'Previously known to affect spermatogenesis', 'Function', 'Previous ref to RNAi working', 'Other species infertility phenotype?'
    ]
    df[cols1] = dfS[cols1]
    df[cols2] = dfS[cols2]
    df[cols3] = dfS[cols3]

    # Export
    df.to_csv('results/core_DM_meiotic_genes_screened.csv')

    print('Done.')