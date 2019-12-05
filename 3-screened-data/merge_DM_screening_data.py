# coding=utf-8
# Author: Rion B Correia
# Date: Jul 17, 2019
#
# Description: Merges DM Selected Gens with DM Screening Data
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists

def map_multiple_ids(x, d):
    x = x.split(',')
    return ','.join([d[i] for i in x])

if __name__ == '__main__':

    # pipeline = 'all3-conserved'  # 'all3-conserved' or 'all3-pooling-DM'

    # Screened Data (From Experimental Analysis)
    dfSc = pd.read_csv('data/conserved_DM_screened_2019-11-22.csv', index_col=0)
    dfSp = pd.read_csv('data/pooling_DM_screened_2019-11-22.csv', index_col=0)
    dfS = pd.concat([dfSc, dfSp], axis='index', join='outer').drop_duplicates()

    # DGE Genes per pipeline
    dfC = pd.read_csv('../2-core_genes/results/all3-conserved/DM_meiotic_genes.csv', index_col=0)
    dfP = pd.read_csv('../2-core_genes/results/all3-pooling-DM/DM_meiotic_genes.csv', index_col=0)

    ####
    # Meta Genes (FOR ANALYSIS SIMPLICITY)
    dfM = pd.read_csv('../2-core_genes/results/all3-pooling-DM/meta_meiotic_genes.csv', index_col='id_eggnog')
    dfM_HS = pd.read_csv('../2-core_genes/results/all3-pooling-DM/HS_meiotic_genes.csv', index_col='id_string', usecols=['id_gene', 'id_string', 'gene'])
    dfM_MM = pd.read_csv('../2-core_genes/results/all3-pooling-DM/MM_meiotic_genes.csv', index_col='id_string', usecols=['id_gene', 'id_string', 'gene'])
    dfM_DM = pd.read_csv('../2-core_genes/results/all3-pooling-DM/DM_meiotic_genes.csv', index_col='id_string', usecols=['id_gene', 'id_string', 'gene'])

    dfM['id_gene_HS'] = dfM['id_string_HS'].apply(map_multiple_ids, args=(dfM_HS['id_gene'].to_dict(),))
    dfM['id_gene_MM'] = dfM['id_string_MM'].apply(map_multiple_ids, args=(dfM_MM['id_gene'].to_dict(),))
    dfM['id_gene_DM'] = dfM['id_string_DM'].apply(map_multiple_ids, args=(dfM_DM['id_gene'].to_dict(),))

    dfM['gene_HS'] = dfM['id_string_HS'].apply(map_multiple_ids, args=(dfM_HS['gene'].to_dict(),))
    dfM['gene_MM'] = dfM['id_string_MM'].apply(map_multiple_ids, args=(dfM_MM['gene'].to_dict(),))
    dfM['gene_DM'] = dfM['id_string_DM'].apply(map_multiple_ids, args=(dfM_DM['gene'].to_dict(),))

    dfM['id_gene_DM'] = dfM['id_gene_DM'].str.split(',')
    dfM = dfM.explode('id_gene_DM').set_index('id_gene_DM') # EXPLODE IS A SUPER COOL FUNCTION!
    ####

    # subset the pooling not to contain the conserved
    dfP = dfP.loc[(~dfP.index.isin(dfC.index)), :]

    cols1 = [
        'Stock Number', 'Transgene', 'status'
    ]
    cols2 = [
        'FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched',
    ]
    cols3 = [
        'Recorded cellular phenotype', 'Previously known to affect male fertility/sperm cells', 'Function', 'Previous ref to RNAi working', 'Other species infertility phenotype?', 'Flybase link'
    ]
    # Extra Information from Mouse & Human
    cols4 = [
        'id_gene_HS', 'id_gene_MM', 'gene_HS', 'gene_MM'
    ]
    for df in [dfC, dfP]:
        df[cols1] = dfS[cols1]
        df[cols2] = dfS[cols2]
        df[cols3] = dfS[cols3]
        # Extra Information from Mouse & Human
        df[cols4] = dfM[cols4]

        ####

        ####

    # Export Core
    wdfCfile = 'results/all3-conserved/DM_meiotic_genes_screened.csv'
    ensurePathExists(wdfCfile)
    dfC.to_csv(wdfCfile)

    # Export Pooling
    wdfPfile = 'results/all3-pooling-DM/DM_meiotic_genes_screened.csv'
    ensurePathExists(wdfPfile)
    dfP.to_csv(wdfPfile)

    print('Done.')
