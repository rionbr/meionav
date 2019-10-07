# coding=utf-8
# Author: Rion B Correia
# Date: Jul 17, 2019
#
# Description: Merges DM Selected Gens with DM Screening Data
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists


if __name__ == '__main__':

    pipeline = 'all3-pooling-DM-FDRp01'  # 'all3-conserved-FDRp05' or 'all3-pooling-DM-FDRp01'

    # This is the 'dfU' object from [3]
    rDMfile = '../2-core_genes/results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline)
    df = pd.read_csv(rDMfile, index_col=0)
    
    # Screened Data (From Experimental Analysis)
    dfS = pd.read_csv('data/core_DM_screened_2019-10-03.csv', index_col=0)

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
    wDMfile = 'results/{pipeline:s}/DM_meiotic_genes_screened.csv'.format(pipeline=pipeline)
    ensurePathExists(wDMfile)
    df.to_csv(wDMfile)

    print('Done.')