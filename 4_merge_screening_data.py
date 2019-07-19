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
    df_up = pd.read_csv('results/core_DM_meiotic_genes_up.csv').set_index('id_gene_DM')
    df_down = pd.read_csv('results/core_DM_meiotic_genes_down.csv').set_index('id_gene_DM')
    
    # Exclusion (genes are unique among the three tables)
    df_up_exc = df_up.loc[ ~df_up.index.isin(df_down.index), :].copy()
    df_down_exc = df_down.loc[ ~df_down.index.isin(df_up.index), :].copy()
    # Intersection between up and down
    df_updown_exc = pd.merge(df_up, df_down, how='inner', left_index=True, right_index=True)

    # Screened Data (From Paulo)
    df_S_up = pd.read_csv('results/screening/DM_MeioticScreen_Up.csv').set_index('id_gene_DM')
    df_S_down = pd.read_csv('results/screening/DM_MeioticScreen_Down.csv').set_index('id_gene_DM')
    df_S_updown = pd.read_csv('results/screening/DM_MeioticScreen_UpDown.csv').set_index('id_gene_DM')


    for ud, df, df_S in zip(['up', 'down', 'updown'], [df_up_exc, df_down_exc, df_updown_exc], [df_S_up, df_S_down, df_S_updown]):
        cols1 = [
            'StockNumber', 'FlyBase Genotype'
        ]
        cols2 = [
            'FT1 eggs', 'FT1 hatched', 'FT2 eggs', 'FT2 hatched', 'FT3 eggs', 'FT3 hatched', 'FT4 eggs', 'FT4 hatched',
        ]
        cols3 = [
            'Phenotype or literature', 'Previously known to affect spermatogenesis', 'Function', 'Overlapping', 'Previous ref to RNAi working'
        ]
        df[cols1] = df_S[cols1]
        df[cols2] = df_S[cols2]
        df[cols3] = df_S[cols3]
        df.to_csv('results/screened/core_DM_screened_{:s}.csv'.format(ud))

    print('Done.')