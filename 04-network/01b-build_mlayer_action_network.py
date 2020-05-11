# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Builds a MultiLayer Action (directed) network (HS, MM & DM) based on genes found by DGE with StringDB edges.
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from itertools import chain, product
from utils import ensurePathExists


def actionapply(dft):
    mode = dft['mode'].tolist()
    action = dft['action'].tolist()
    score = dft['score'].mean()
    return pd.Series({
        # mode
        'reaction': True if 'reaction' in mode else np.nan,
        'expression': True if 'expression' in mode else np.nan,
        'activation': True if 'activation' in mode else np.nan,
        'ptmod': True if 'ptmod' in mode else np.nan,
        'binding': True if 'binding' in mode else np.nan,
        'catalysis': True if 'catalysis' in mode else np.nan,
        # action
        'inhibition': True if 'inhibition' in action else np.nan,
        'activation': True if 'activation' in action else np.nan,
        #
        'weight': score / 1000
    })


if __name__ == '__main__':

    G = nx.DiGraph()
    minLogTPM = math.log2(2)


    ##
    # HS Network Data
    ##
    print('Processing HS data')
    df_HS = pd.read_csv("../2-core_genes/results/HS-FPKM_genes.csv.gz", index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'FPKM', 'TPM', 'Spermatocytes', 'biotype'])
    df_HS = df_HS.rename(columns={'gene': 'label'})  # Rename
    # Only TPM > log(2) & Spermatocyes
    df_HS = df_HS.loc[((df_HS['TPM'] >= minLogTPM) & (df_HS['Spermatocytes'] == True)), :]
    # Identify DGE in Meiotic Entry/Exit Genes
    df_HS_DGE = pd.read_csv("../2-core_genes/results/HS-DE_genes.csv.gz", index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'Cyte_vs_Gonia', 'Tid_vs_Cyte'])
    df_HS['meiotic-entry'] = df_HS_DGE['Cyte_vs_Gonia']
    df_HS['meiotic-exit'] = df_HS_DGE['Tid_vs_Cyte']
    # Map id_gene -> id_string
    dict_HS_id_gene_to_id_string = df_HS.explode('id_string').reset_index().set_index('id_string')['id_gene'].to_dict()
    # Add Core and Mammal information
    s_HS_core = pd.read_csv("../2-core_genes/results/pipeline-core/HS_meiotic_genes.csv", index_col='id_gene', usecols=['id_gene', 'id_string'])
    s_HS_mammals = pd.read_csv("../2-core_genes/results/pipeline-mammals/HS_meiotic_genes.csv", index_col='id_gene', usecols=['id_gene', 'id_string'])
    # Identify Core
    df_HS['core'] = df_HS.index.map(lambda x: True if x in s_HS_core.index else np.nan)
    df_HS['mammals'] = df_HS.index.map(lambda x: True if x in s_HS_mammals.index else np.nan)
    # Bag of String (some genes have >1 id_string)
    set_HS_id_strings = set(np.hstack(df_HS['id_string'].dropna()))
    # logFPKM
    df_HS['logFPKM'] = df_HS['FPKM'].apply(lambda x: np.log2(x + 1))
    df_HS.drop(['FPKM'], axis='columns', inplace=True)
    # Identify layer
    df_HS['layer'] = 'HS'

    # Add tested results
    print('Adding laboratory results')
    df_HS_S = pd.read_csv('../3-screened-data/data/literature_HS_2020-01-07.csv', index_col='id_gene_HS')
    df_HS['known-HS-phenotype'] = df_HS_S['Others HS pheno code']

    # Debug
    # df_HS = df_HS.loc[(df_HS['core'] == True), :]
    # print(df_HS.head())

    # Add Nodes with tuple(id, attrs)
    idxs = df_HS.index.to_list()
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_HS.to_dict(orient='rows')]
    G.add_nodes_from([(i, d) for i, d in zip(idxs, attrs)])

    #
    # HS Actions
    #
    print('Adding actions')
    df_HS_actions = pd.read_csv("../StringDB/9606/9606.protein.actions.v11.0.txt.gz", sep='\t')
    # Rename id_string columns
    df_HS_actions = df_HS_actions.rename(columns={'item_id_a': 'id_string_i', 'item_id_b': 'id_string_j'})
    # Reduce Search Space
    df_HS_actions = df_HS_actions.loc[(df_HS_actions['id_string_i'].isin(set_HS_id_strings) & df_HS_actions['id_string_j'].isin(set_HS_id_strings)), :]
    # Only considers directional evidence where a_is_acting
    df_HS_actions['is_directional'] = df_HS_actions['is_directional'].map(lambda x: True if x == 't' else False)
    df_HS_actions['a_is_acting'] = df_HS_actions['a_is_acting'].map(lambda x: True if x == 't' else False)
    df_HS_actions = df_HS_actions.loc[(df_HS_actions['is_directional'] == True) & (df_HS_actions['a_is_acting'] == True), :]
    # Group-apply to deal with duplicates
    df_HS_actions = df_HS_actions.groupby(['id_string_i', 'id_string_j']).apply(actionapply)
    df_HS_actions = df_HS_actions.reset_index()
    # Map id_string to id_gene
    df_HS_actions['id_gene_i'] = df_HS_actions['id_string_i'].map(lambda x: dict_HS_id_gene_to_id_string[x])
    df_HS_actions['id_gene_j'] = df_HS_actions['id_string_j'].map(lambda x: dict_HS_id_gene_to_id_string[x])
    # Set Index
    df_HS_actions = df_HS_actions.set_index(['id_gene_i', 'id_gene_j'])
    # Add/Update Edge with tuple(id, attrs)
    idxs = df_HS_actions.index.to_list()
    # edge data (removing NaN values)
    # attrs = df_HS_actions.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_HS_actions.to_dict(orient='rows')]
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])


    ##
    # Mouse Musculus (MM) Network Data
    ##
    print('Processing MM data')
    df_MM = pd.read_csv("../2-core_genes/results/MM-FPKM_genes.csv.gz", index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'FPKM', 'TPM', 'Spermatocytes', 'biotype'])
    df_MM = df_MM.rename(columns={'gene': 'label'})  # Rename
    # Only TPM >= log(2) & Spermatocyes
    df_MM = df_MM.loc[((df_MM['TPM'] >= minLogTPM) & (df_MM['Spermatocytes'] == True)), :]
    # Identify DGE in Meiotic Entry/Exit Genes
    df_MM_DGE = pd.read_csv("../2-core_genes/results/MM-DE_genes.csv.gz", index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'Cyte_vs_Gonia', 'Tid_vs_Cyte'])
    df_MM['meiotic-entry'] = df_MM_DGE['Cyte_vs_Gonia']
    df_MM['meiotic-exit'] = df_MM_DGE['Tid_vs_Cyte']
    # Map id_gene -> id_string
    dict_MM_id_gene_to_id_string = df_MM.explode('id_string').reset_index().set_index('id_string')['id_gene'].to_dict()
    # Add Core and Mammal information#
    s_MM_core = pd.read_csv("../2-core_genes/results/pipeline-core/MM_meiotic_genes.csv", index_col='id_gene', usecols=['id_gene', 'id_string'])
    s_MM_mammals = pd.read_csv("../2-core_genes/results/pipeline-mammals/MM_meiotic_genes.csv", index_col='id_gene', usecols=['id_gene', 'id_string'])
    # Identify Core
    df_MM['core'] = df_MM.index.map(lambda x: True if x in s_MM_core.index else np.nan)
    df_MM['mammals'] = df_MM.index.map(lambda x: True if x in s_MM_mammals.index else np.nan)
    # Bag of String (some genes have >1 id_string)
    set_MM_id_strings = set(np.hstack(df_MM['id_string'].dropna()))
    # logFPKM
    df_MM['logFPKM'] = df_MM['FPKM'].apply(lambda x: np.log2(x + 1))
    df_MM.drop(['FPKM'], axis='columns', inplace=True)
    # Identify layer
    df_MM['layer'] = 'MM'

    # Add tested results
    print('Adding laboratory results')
    df_MM_S = pd.read_csv('../3-screened-data/data/literature_MM_2020-01-07.csv', index_col='id_gene_MM')
    df_MM['known-MM-phenotype'] = df_MM_S['Others MM pheno code']

    # Debug
    # df_MM = df_MM.loc[(df_MM['core'] == True), :]

    # Add Nodes with tuple(id, attrs)
    idxs = df_MM.index.to_list()
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_MM.to_dict(orient='rows')]
    G.add_nodes_from([(i, d) for i, d in zip(idxs, attrs)])


    ##
    # MM Actions
    ##
    df_MM_actions = pd.read_csv("../StringDB/10090/10090.protein.actions.v11.0.txt.gz", sep='\t')
    # Rename id_string columns
    df_MM_actions = df_MM_actions.rename(columns={'item_id_a': 'id_string_i', 'item_id_b': 'id_string_j'})
    # Reduce Search Space
    df_MM_actions = df_MM_actions.loc[(df_MM_actions['id_string_i'].isin(set_MM_id_strings) & df_MM_actions['id_string_j'].isin(set_MM_id_strings)), :]
    # Only considers directional evidence where a_is_acting
    df_MM_actions['is_directional'] = df_MM_actions['is_directional'].map(lambda x: True if x == 't' else False)
    df_MM_actions['a_is_acting'] = df_MM_actions['a_is_acting'].map(lambda x: True if x == 't' else False)
    df_MM_actions = df_MM_actions.loc[(df_MM_actions['is_directional'] == True) & (df_MM_actions['a_is_acting'] == True), :]
    # Group-apply to deal with duplicates
    df_MM_actions = df_MM_actions.groupby(['id_string_i', 'id_string_j']).apply(actionapply)
    df_MM_actions = df_MM_actions.reset_index()
    # Map id_string to id_gene
    df_MM_actions['id_gene_i'] = df_MM_actions['id_string_i'].map(lambda x: dict_MM_id_gene_to_id_string[x])
    df_MM_actions['id_gene_j'] = df_MM_actions['id_string_j'].map(lambda x: dict_MM_id_gene_to_id_string[x])
    # Set Index
    df_MM_actions = df_MM_actions.set_index(['id_gene_i', 'id_gene_j'])
    # Add/Update Edge with tuple(id, attrs)
    idxs = df_MM_actions.index.to_list()
    # edge data (removing NaN values)
    # attrs = df_MM_actions.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list() # THIS IS TOO SLOW
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_MM_actions.to_dict(orient='rows')]
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])


    ##
    # Drosophila Melanogaster (DM) Network Data
    ##
    print('Processing DM data')
    df_DM = pd.read_csv("../2-core_genes/results/DM-FPKM_genes.csv.gz", index_col='id_gene', usecols=['id_gene', 'id_string', 'gene', 'FPKM', 'TPM', 'Middle', 'biotype'])
    df_DM = df_DM.rename(columns={'gene': 'label'})  # Rename
    # Only TPM >= log(2)
    df_DM = df_DM.loc[((df_DM['TPM'] >= minLogTPM) & (df_DM['Middle'] == True)), :]
    # Identify DGE in Meiotic Entry/Exit Genes
    df_DM_DGE = pd.read_csv("../2-core_genes/results/DM-DE_genes.csv.gz", index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'Middle_vs_Apical', 'Basal_vs_Middle'])
    df_DM['meiotic-entry'] = df_DM_DGE['Middle_vs_Apical']
    df_DM['meiotic-exit'] = df_DM_DGE['Basal_vs_Middle']
    # Map id_gene -> id_string
    dict_DM_id_gene_to_id_string = df_DM.explode('id_string').reset_index().set_index('id_string')['id_gene'].to_dict()
    # Add Core information
    s_DM_core = pd.read_csv("../2-core_genes/results/pipeline-core/DM_meiotic_genes.csv", index_col='id_gene', usecols=['id_gene', 'id_string'])
    # Identify Core
    df_DM['core'] = df_DM.index.map(lambda x: True if x in s_DM_core.index else np.nan)
    # Bag of String (some genes have >1 id_string)
    set_DM_id_strings = set(np.hstack(df_DM['id_string'].dropna()))
    # logFPKM
    df_DM['logFPKM'] = df_DM['FPKM'].apply(lambda x: np.log2(x + 1))
    df_DM.drop(['FPKM'], axis='columns', inplace=True)
    # Identify layer
    df_DM['layer'] = 'DM'

    # Add tested results
    print('Adding laboratory results')
    df_DM_S = pd.read_csv('../3-screened-data/data/core_DM_screened_2020-02-14.csv', index_col='id_gene')
    # Calculations
    for ft in range(1, 5):
        col_eggs = 'FT{:d} eggs'.format(ft)
        col_hatched = 'FT{:d} hatched'.format(ft)
        col_fertate = 'FT{:d} fert-rate'.format(ft)
        df_DM_S[col_fertate] = df_DM_S[col_hatched] / df_DM_S[col_eggs]
    df_DM['mean-fert-rate'] = df_DM_S[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].mean(axis=1)
    df_DM['std-fert-rate'] = df_DM_S[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].std(axis=1)
    df_DM[['status', 'recorded-phenotype', 'previously-reported', 'function', 'validated-rnai', 'new-DM-phenotype', 'known-DM-phenotype']] = df_DM_S[['Status', 'Recorded cellular phenotype', 'Previously reported in DM?', 'Function', 'Previous ref to RNAi working?', 'Our DM pheno code', 'Others DM pheno code']]
    # Debug
    # df_DM = df_DM.loc[(df_DM['core'] == True), :]
    # print(df_DM.head())

    # Add Nodes with tuple(id, attrs)
    idxs = df_DM.index.to_list()
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_DM.to_dict(orient='rows')]
    G.add_nodes_from([(i, d) for i, d in zip(idxs, attrs)])

    ##
    # DM Actions
    ##
    df_DM_actions = pd.read_csv("../StringDB/7227/7227.protein.actions.v11.0.txt.gz", sep='\t')
    # Rename id_string columns
    df_DM_actions = df_DM_actions.rename(columns={'item_id_a': 'id_string_i', 'item_id_b': 'id_string_j'})
    # Reduce Search Space
    df_DM_actions = df_DM_actions.loc[(df_DM_actions['id_string_i'].isin(set_DM_id_strings) & df_DM_actions['id_string_j'].isin(set_DM_id_strings)), :]
    # Only considers directional evidence where a_is_acting
    df_DM_actions['is_directional'] = df_DM_actions['is_directional'].map(lambda x: True if x == 't' else False)
    df_DM_actions['a_is_acting'] = df_DM_actions['a_is_acting'].map(lambda x: True if x == 't' else False)
    df_DM_actions = df_DM_actions.loc[(df_DM_actions['is_directional'] == True) & (df_DM_actions['a_is_acting'] == True), :]
    # Groupby mode
    df_DM_actions = df_DM_actions.groupby(['id_string_i', 'id_string_j']).apply(actionapply)
    df_DM_actions = df_DM_actions.reset_index()
    # Map id_string to id_gene
    df_DM_actions['id_gene_i'] = df_DM_actions['id_string_i'].map(lambda x: dict_DM_id_gene_to_id_string[x])
    df_DM_actions['id_gene_j'] = df_DM_actions['id_string_j'].map(lambda x: dict_DM_id_gene_to_id_string[x])
    # Set Index
    df_DM_actions = df_DM_actions.set_index(['id_gene_i', 'id_gene_j'])
    # Add/Update Edge with tuple(id, attrs)
    idxs = df_DM_actions.index.to_list()
    # edge data (removing NaN values)
    # attrs = df_DM_actions.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()  # THIS IS TOO SLOW
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_DM_actions.to_dict(orient='rows')]
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])
