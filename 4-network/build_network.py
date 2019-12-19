# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Builds a MultiLayer network (HS, MM & DM) based on genes found by DGE with StringDB edges.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from itertools import chain, product
from utils import ensurePathExists
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors


cmap_meanfertrate = colors.LinearSegmentedColormap.from_list(name='cmap-mean-fert-rate', colors=['#d62728', '#1f77b4'], N=256)


def fert_rate_color(x):
    if pd.isnull(x):
        return '#FFFFFF'  # white
    else:
        return colors.to_hex(cmap_meanfertrate(x))  # color


if __name__ == '__main__':

    G = nx.Graph()

    ##
    # HS Network Data
    ##
    """
    print('Processing HS data')
    df_HS = pd.read_csv("../2-core_genes/results/HS-DE_genes.csv", index_col='id_string', usecols=['id_string', 'id_gene', 'gene', 'Cyte_vs_Gonia', 'Tid_vs_Cyte', 'biotype'])
    df_HS = df_HS.rename(columns={'gene':'label'})  # Rename
    # Add Core and Pool information#
    s_HS_core = pd.read_csv("../2-core_genes/results/all3-conserved/HS_meiotic_genes.csv", index_col='id_string', usecols=['id_string'])
    s_HS_pool = pd.read_csv("../2-core_genes/results/all3-pooling-DM/HS_meiotic_genes.csv", index_col='id_string', usecols=['id_string'])
    #
    df_HS['conserved'] = df_HS.index.map(lambda x: True if x in s_HS_core.index else np.nan)
    df_HS['pooling'] = df_HS.index.map(lambda x: True if x in s_HS_pool.index else np.nan)
    # Identify layer
    df_HS['layer'] = 'HS'

    # Debug
    df_HS = df_HS.loc[df_HS['conserved'] == True, :]
    print(df_HS.head())

    # Add Nodes with tuple(id, attrs)
    G.add_nodes_from(zip(df_HS.index, df_HS.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()))

    #
    # HS Links
    #
    df_HS_links = pd.read_csv("../StringDB/9606/9606.protein.links.full.v11.0.txt.gz", sep=' ', index_col=[0, 1])
    # Reduce Search Space
    df_HS_links = df_HS_links.loc[(df_HS_links.index.isin(df_HS.index, level=0) & df_HS_links.index.isin(df_HS.index, level=1)), :]
    # replace 0 for NaN
    df_HS_links = df_HS_links.replace({0.0: np.nan})
    # Add edge type
    df_HS_links['type'] = 'intra'

    # Add Edges with tuple(id, attrs)
    # edge i,j indexes
    idxs = df_HS_links.index.to_list()
    # edge data (removing NaN values)
    attrs = df_HS_links.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])

    #
    # HS Actions
    #
    df_HS_actions = pd.read_csv("../StringDB/9606/9606.protein.actions.v11.0.txt.gz", sep='\t', index_col=[0, 1])
    # Reduce Search Space
    df_HS_actions = df_HS_actions.loc[(df_HS_actions.index.isin(df_HS.index, level=0) & df_HS_actions.index.isin(df_HS.index, level=1)), :]
    # replace 0 for NaN
    df_HS_actions = df_HS_actions.replace({0.0: np.nan})

    # Add/Update Edge with tuple(id, attrs)
    idx = df_HS_actions.index.to_list()
    # edge data (removing NaN values)
    attrs = df_HS_actions.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])

    ##
    # Mouse Musculus (MM) Network Data
    ##
    print('Processing MM data')
    df_MM = pd.read_csv("../2-core_genes/results/MM-DE_genes.csv", index_col='id_string', usecols=['id_string', 'id_gene', 'gene', 'Cyte_vs_Gonia', 'Tid_vs_Cyte', 'biotype'])
    df_MM = df_MM.rename(columns={'gene':'label'})  # Rename
    # Add Core and Pool information#
    s_MM_core = pd.read_csv("../2-core_genes/results/all3-conserved/MM_meiotic_genes.csv", index_col='id_string', usecols=['id_string'])
    s_MM_pool = pd.read_csv("../2-core_genes/results/all3-pooling-DM/MM_meiotic_genes.csv", index_col='id_string', usecols=['id_string'])
    #
    df_MM['conserved'] = df_MM.index.map(lambda x: True if x in s_MM_core.index else np.nan)
    df_MM['pooling'] = df_MM.index.map(lambda x: True if x in s_MM_pool.index else np.nan)
    # Identify layer
    df_MM['layer'] = 'MM'

    # Debug
    df_MM = df_MM.loc[df_MM['conserved'] == True, :]

    # Add Nodes with tuple(id, attrs)
    G.add_nodes_from(zip(df_MM.index, df_MM.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()))

    ##
    # MM Links
    ##
    df_MM_links = pd.read_csv("../StringDB/10090/10090.protein.links.full.v11.0.txt.gz", sep=' ', index_col=[0, 1])
    # Reduce Search Space
    df_MM_links = df_MM_links.loc[(df_MM_links.index.isin(df_MM.index, level=0) & df_MM_links.index.isin(df_MM.index, level=1)), :]
    # replace 0 for NaN
    df_MM_links = df_MM_links.replace({0.0: np.nan})
    # Add edge type
    df_MM_links['type'] = 'intra'
    

    # Add Edges with tuple(id, attrs)
    # edge i,j indexes
    idxs = df_MM_links.index.to_list()
    # edge data (removing NaN values)
    attrs = df_MM_links.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])
    
    ##
    # MM Actions
    ##
    df_MM_actions = pd.read_csv("../StringDB/10090/10090.protein.actions.v11.0.txt.gz", sep='\t', index_col=[0, 1])
    # Reduce Search Space
    df_MM_actions = df_MM_actions.loc[(df_MM_actions.index.isin(df_MM.index, level=0) & df_MM_actions.index.isin(df_MM.index, level=1)), :]
    # replace 0 for NaN
    df_MM_actions = df_MM_actions.replace({0.0: np.nan})

    # Add/Update Edge with tuple(id, attrs)
    idx = df_MM_actions.index.to_list()
    # edge data (removing NaN values)
    attrs = df_MM_actions.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])
    """
    ##
    # Drosophila Melanogaster (DM) Network Data
    ##
    print('Processing DM data')
    df_DM = pd.read_csv("../2-core_genes/results/DM-DE_genes.csv", index_col='id_string', usecols=['id_string', 'id_gene', 'gene', 'Middle_vs_Apical', 'Basal_vs_Middle', 'biotype'])
    df_DM = df_DM.rename(columns={'gene':'label'})  # Rename
    # Add Core and Pool information#
    s_DM_core = pd.read_csv("../2-core_genes/results/all3-conserved/DM_meiotic_genes.csv", index_col='id_string', usecols=['id_string'])
    s_DM_pool = pd.read_csv("../2-core_genes/results/all3-pooling-DM/DM_meiotic_genes.csv", index_col='id_string', usecols=['id_string'])
    #
    df_DM['conserved'] = df_DM.index.map(lambda x: True if x in s_DM_core.index else False)
    df_DM['pooling'] = df_DM.index.map(lambda x: True if x in s_DM_pool.index else False)

    # Add tested results
    print('Adding laboratory results')
    df_DM_Sc = pd.read_csv('../3-screened-data/data/conserved_DM_screened_2019-11-22.csv', index_col='id_gene')
    df_DM_Sp = pd.read_csv('../3-screened-data/data/pooling_DM_screened_2019-11-22.csv', index_col='id_gene')
    df_DM_S = pd.concat([df_DM_Sc, df_DM_Sp], axis='index', join='outer').drop_duplicates()
    # Map id_string
    DMmap = df_DM.reset_index().set_index('id_gene')['id_string'].to_dict()
    df_DM_S['id_string'] = df_DM_S.index.map(DMmap)
    df_DM_S = df_DM_S.reset_index().set_index('id_string')
    # Calculations
    for ft in range(1, 5):
        col_eggs = 'FT{:d} eggs'.format(ft)
        col_hatched = 'FT{:d} hatched'.format(ft)
        col_fertate = 'FT{:d} fert-rate'.format(ft)
        df_DM_S[col_fertate] = df_DM_S[col_hatched] / df_DM_S[col_eggs]
    df_DM['mean fert-rate'] = df_DM_S[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].mean(axis=1)
    df_DM['std fert-rate'] = df_DM_S[['FT1 fert-rate', 'FT2 fert-rate', 'FT3 fert-rate', 'FT4 fert-rate']].std(axis=1)
    # Identify layer
    df_DM['layer'] = 'DM'
    # Node Color (baed on Mean Fert-Rate)
    df_DM['color'] = df_DM['mean fert-rate'].apply(fert_rate_color)
    # Size
    df_DM['size'] = 25

    # Debug
    df_DM = df_DM.loc[df_DM['conserved'] == True, :]
    print(df_DM.head())

    # Add Nodes with tuple(id, attrs)
    G.add_nodes_from(zip(df_DM.index, df_DM.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()))

    ##
    # DS Links
    ##
    df_DM_links = pd.read_csv("../StringDB/7227/7227.protein.links.full.v11.0.txt.gz", sep=' ', index_col=[0, 1])
    # Reduce Search Space
    df_DM_links = df_DM_links.loc[(df_DM_links.index.isin(df_DM.index, level=0) & df_DM_links.index.isin(df_DM.index, level=1)), :]
    # replace 0 for NaN
    df_DM_links = df_DM_links.replace({0.0: np.nan})
    # Add edge type
    df_DM_links['type'] = 'intra'
    # Weight
    df_DM_links['weight'] = df_DM_links['combined_score'] / 1000

    # Debug
    print(df_DM_links.head())

    # Add Edges with tuple(id, attrs)
    # edge i,j indexes
    idxs = df_DM_links.index.to_list()
    # edge data (removing NaN values)
    attrs = df_DM_links.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])

    ##
    # DM Actions
    ##
    df_DM_actions = pd.read_csv("../StringDB/7227/7227.protein.actions.v11.0.txt.gz", sep='\t', index_col=[0, 1])
    # Reduce Search Space
    df_DM_actions = df_DM_actions.loc[(df_DM_actions.index.isin(df_DM.index, level=0) & df_DM_actions.index.isin(df_DM.index, level=1)), :]
    # replace 0 for NaN
    df_DM_actions = df_DM_actions.replace({0.0: np.nan})

    # Debug
    print(df_DM_actions.head())

    # Add/Update Edge with tuple(id, attrs)
    idx = df_DM_actions.index.to_list()
    # edge data (removing NaN values)
    attrs = df_DM_actions.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])

    ##
    # Add cross-layer edges
    ##
    """
    dfM = pd.read_csv('../2-core_genes/results/meta_meiotic_genes.csv', index_col='id_eggnog')
    # Change column type from string to list
    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    cross_edges = []
    for id_meta, row in dfM.iterrows():
        id_string_HS = row['id_string_HS']
        id_string_MM = row['id_string_MM']
        id_string_DM = row['id_string_DM']
        # all pairs for each pairwise product
        all_pairs = chain(product(*[id_string_HS, id_string_MM]), product(*[id_string_HS, id_string_DM]), product(*[id_string_MM, id_string_DM]))
        cross_edges.extend(all_pairs)
    G.add_edges_from(cross_edges, type='cross')
    """
    ##
    # Export
    ##
    wGfile_gpickle = 'results/complete_mlayer_graph.gpickle'
    ensurePathExists(wGfile_gpickle)
    nx.write_gpickle(G, wGfile_gpickle)

    wGfile_graphml = 'results/complete_mlayer_graph.graphml'
    nx.write_graphml(G, wGfile_graphml)
