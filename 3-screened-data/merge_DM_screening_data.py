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

    # pipeline = 'all3-conserved'  # 'all3-conserved' or 'all3-pooling-DM'

    # Screened Data (From Experimental Analysis)
    dfSc = pd.read_csv('data/conserved_DM_screened_2019-11-22.csv', index_col=0)
    dfSp = pd.read_csv('data/pooling_DM_screened_2019-11-22.csv', index_col=0)
    dfS = pd.concat([dfSc, dfSp], axis='index', join='outer').drop_duplicates()

    # DGE Genes per pipeline
    dfC = pd.read_csv('../2-core_genes/results/all3-conserved/DM_meiotic_genes.csv', index_col=0)
    dfP = pd.read_csv('../2-core_genes/results/all3-pooling-DM/DM_meiotic_genes.csv', index_col=0)

    # subset the pooling not to contain the conserved
    dfP = dfP.loc[ (~dfP.index.isin(dfC.index)), :]

    cols1 = [
        'Stock Number', 'Transgene', 'status'
    ]
    cols2 = [
        'FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched',
    ]
    cols3 = [
        'Recorded cellular phenotype', 'Previously known to affect male fertility/sperm cells', 'Function', 'Previous ref to RNAi working', 'Other species infertility phenotype?', 'Flybase link'
    ]
    for df in [dfC, dfP]:
        df[cols1] = dfS[cols1]
        df[cols2] = dfS[cols2]
        df[cols3] = dfS[cols3]

    # Export Core
    wdfCfile = 'results/all3-conserved/DM_meiotic_genes_screened.csv'
    ensurePathExists(wdfCfile)
    dfC.to_csv(wdfCfile)

    # Export Pooling
    wdfPfile = 'results/all3-pooling-DM/DM_meiotic_genes_screened.csv'
    ensurePathExists(wdfPfile)
    dfP.to_csv(wdfPfile)

    print('Done.')