# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Builds a MultiLayer network (HS, MM & DM) based on genes found by DGE with StringDB edges.
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
from collections import defaultdict, Counter
import argparse


def transpose_variable_across_layers(G, variable, combination='sum'):
    """ Method to transpose results in one layer to another. Note duplicate results are either summed or set to majority."""
    dict_i_values = {i: d[variable] for i, d in G.nodes(data=True) if d.get(variable, None) is not None}
    dict_j_values = defaultdict(list)
    for i, v in dict_i_values.items():
        cross_edges = [j for _i, j, d in G.edges(i, data=True) if d.get('type', None) == 'cross']
        for j in cross_edges:
            dict_j_values[j].append(v)
    # Combine multiple values
    if combination == 'sum':
        dict_j_values = {k: sum(l) for k, l in dict_j_values.items()}
    elif combination == 'majority':
        dict_j_values = {k: Counter(l).most_common()[0][0] for k, l in dict_j_values.items()}
    else:
        TypeError("Combination must be either 'sum', or 'majority'.")
    # Set attributes to network
    nx.set_node_attributes(G, values=dict_j_values, name=variable)
    return G


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--minLogTPM", default=2, type=int, help="minLogTPM = math.log2(x). Defaults to 2.")
    args = parser.parse_args()
    #
    minLogTPM = math.log2(args.minLogTPM)
    celltype = args.celltype  # spermatocyte or enterocyte

    #
    # Init
    #
    G = nx.Graph()

    ##
    # HS Network Data
    ##
    print('Processing HS data')
    df_HS = pd.read_csv("../2-core_genes/results/HS-FPKM-{celltype:s}.csv.gz".format(celltype=celltype), index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'FPKM', 'TPM', celltype, 'biotype'])
    df_HS = df_HS.rename(columns={'gene': 'label'})  # Rename
    # Only TPM > log(2) & Spermatocyes
    df_HS = df_HS.loc[((df_HS['TPM'] >= minLogTPM) & (df_HS[celltype] == True)), :]
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
    # HS Links
    #
    print('Adding links')
    df_HS_links = pd.read_csv("../data/StringDB/9606/9606.protein.links.full.v11.0.txt.gz", sep=' ', usecols=['protein1', 'protein2', 'experiments', 'combined_score'])
    # Rename id_string columns
    df_HS_links = df_HS_links.rename(columns={'protein1': 'id_string_i', 'protein2': 'id_string_j'})
    # Reduce Search Space
    df_HS_links = df_HS_links.loc[(df_HS_links['id_string_i'].isin(set_HS_id_strings) & df_HS_links['id_string_j'].isin(set_HS_id_strings)), :]
    # Map id_string to id_gene
    df_HS_links['id_gene_i'] = df_HS_links['id_string_i'].map(lambda x: dict_HS_id_gene_to_id_string[x])
    df_HS_links['id_gene_j'] = df_HS_links['id_string_j'].map(lambda x: dict_HS_id_gene_to_id_string[x])
    # replace 0 for NaN
    df_HS_links = df_HS_links.replace({0.0: np.nan})
    # Add edge type
    df_HS_links['type'] = 'intra'
    # Weight
    df_HS_links['weight'] = df_HS_links['combined_score'] / 1000
    # Set Index
    df_HS_links = df_HS_links.set_index(['id_gene_i', 'id_gene_j'])
    # Add/Updates Edges with tuple(id, attrs)
    idxs = df_HS_links.index.to_list()
    # edge data (removing NaN values)
    # attrs = df_HS_links.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()  # THIS IS TOO SLOW
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_HS_links.to_dict(orient='rows')]  # THIS IS FASTER, BUT STILL SLOW
    # attrs = df_HS_links.to_dict(orient='rows')
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])

    ##
    # Mouse Musculus (MM) Network Data
    ##
    print('Processing MM data')
    df_MM = pd.read_csv("../2-core_genes/results/MM-FPKM-{celltype:s}.csv.gz".format(celltype=celltype), index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'FPKM', 'TPM', celltype, 'biotype'])
    df_MM = df_MM.rename(columns={'gene': 'label'})  # Rename
    # Only TPM >= log(2) & Spermatocyes
    df_MM = df_MM.loc[((df_MM['TPM'] >= minLogTPM) & (df_MM[celltype] == True)), :]
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
    # MM Links
    ##
    print('Adding links')
    df_MM_links = pd.read_csv("../data/StringDB/10090/10090.protein.links.full.v11.0.txt.gz", sep=' ', usecols=['protein1', 'protein2', 'experiments', 'combined_score'])
    # Rename id_string columns
    df_MM_links = df_MM_links.rename(columns={'protein1': 'id_string_i', 'protein2': 'id_string_j'})
    # Reduce Search Space
    df_MM_links = df_MM_links.loc[(df_MM_links['id_string_i'].isin(set_MM_id_strings) & df_MM_links['id_string_j'].isin(set_MM_id_strings)), :]
    # Map id_string to id_gene
    df_MM_links['id_gene_i'] = df_MM_links['id_string_i'].map(lambda x: dict_MM_id_gene_to_id_string[x])
    df_MM_links['id_gene_j'] = df_MM_links['id_string_j'].map(lambda x: dict_MM_id_gene_to_id_string[x])
    # replace 0 for NaN
    df_MM_links = df_MM_links.replace({0.0: np.nan})
    # Add edge type
    df_MM_links['type'] = 'intra'
    # Weight
    df_MM_links['weight'] = df_MM_links['combined_score'] / 1000
    # Set Index
    df_MM_links = df_MM_links.set_index(['id_gene_i', 'id_gene_j'])
    # Add/Update Edges with tuple(id, attrs)
    idxs = df_MM_links.index.to_list()
    # edge data (removing NaN values)
    # attrs = df_MM_links.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()  # THIS IS TOO SLOW
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_MM_links.to_dict(orient='rows')]
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])

    ##
    # Drosophila Melanogaster (DM) Network Data
    ##
    print('Processing DM data')
    df_DM = pd.read_csv("../2-core_genes/results/DM-FPKM-{celltype:s}.csv.gz".format(celltype=celltype), index_col='id_gene', usecols=['id_gene', 'id_string', 'gene', 'FPKM', 'TPM', celltype, 'biotype'])
    df_DM = df_DM.rename(columns={'gene': 'label'})  # Rename
    # Only TPM >= log(2)
    df_DM = df_DM.loc[((df_DM['TPM'] >= minLogTPM) & (df_DM[celltype] == True)), :]
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
    # DS Links
    ##
    print('Adding links')
    df_DM_links = pd.read_csv("../data/StringDB/7227/7227.protein.links.full.v11.0.txt.gz", sep=' ', usecols=['protein1', 'protein2', 'experiments', 'combined_score'])
    # Rename id_string columns
    df_DM_links = df_DM_links.rename(columns={'protein1': 'id_string_i', 'protein2': 'id_string_j'})
    # Reduce Search Space
    df_DM_links = df_DM_links.loc[(df_DM_links['id_string_i'].isin(set_DM_id_strings) & df_DM_links['id_string_j'].isin(set_DM_id_strings)), :]
    # Map id_string to id_gene
    df_DM_links['id_gene_i'] = df_DM_links['id_string_i'].map(lambda x: dict_DM_id_gene_to_id_string[x])
    df_DM_links['id_gene_j'] = df_DM_links['id_string_j'].map(lambda x: dict_DM_id_gene_to_id_string[x])
    # replace 0 for NaN
    df_DM_links = df_DM_links.replace({0.0: np.nan})
    # Add edge type
    df_DM_links['type'] = 'intra'
    # Weight
    df_DM_links['weight'] = df_DM_links['combined_score'] / 1000
    # Set Index
    df_DM_links = df_DM_links.set_index(['id_gene_i', 'id_gene_j'])
    # Add/Update Edges with tuple(id, attrs)
    idxs = df_DM_links.index.to_list()
    # edge data (removing NaN values)
    # attrs = df_DM_links.apply(lambda x: x.dropna().to_dict(), axis='columns').to_list()  # THIS IS TOO SLOW
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_DM_links.to_dict(orient='rows')]
    # Add edges
    G.add_edges_from([(i, j, d) for (i, j), d in zip(idxs, attrs)])

    ##
    # Add cross-layer edges
    ##
    print('Adding cross-layer edges')
    dfM = pd.read_csv('../2-core_genes/results/meta_meiotic_genes.csv', index_col='id_eggnog')
    # Change column type from string to list
    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    set_all_node_ids = set(G.nodes())
    cross_edges = []
    for id_meta, row in dfM.iterrows():
        id_string_HS = row['id_string_HS']
        id_string_MM = row['id_string_MM']
        id_string_DM = row['id_string_DM']
        # Map id_string to id_gene
        id_gene_HS = [dict_HS_id_gene_to_id_string[n] for n in id_string_HS if n in dict_HS_id_gene_to_id_string]
        id_gene_MM = [dict_MM_id_gene_to_id_string[n] for n in id_string_MM if n in dict_MM_id_gene_to_id_string]
        id_gene_DM = [dict_DM_id_gene_to_id_string[n] for n in id_string_DM if n in dict_DM_id_gene_to_id_string]
        # only ids already in graph
        id_gene_HS = [n for n in id_gene_HS if n in set_all_node_ids]
        id_gene_MM = [n for n in id_gene_MM if n in set_all_node_ids]
        id_gene_DM = [n for n in id_gene_DM if n in set_all_node_ids]
        # all pairs for each pairwise product
        all_pairs = chain(product(*[id_gene_HS, id_gene_MM]), product(*[id_gene_HS, id_gene_DM]), product(*[id_gene_MM, id_gene_DM]))
        cross_edges.extend(all_pairs)
    G.add_edges_from(cross_edges, type='cross')

    ##
    # Transposing DM Mean Fertility Rate to cross layers
    ##
    print('Tranposing variable across layers')
    G = transpose_variable_across_layers(G, 'mean-fert-rate', combination='sum')
    G = transpose_variable_across_layers(G, 'new-DM-phenotype', combination='majority')
    G = transpose_variable_across_layers(G, 'known-HS-phenotype', combination='majority')
    G = transpose_variable_across_layers(G, 'known-MM-phenotype', combination='majority')
    G = transpose_variable_across_layers(G, 'known-DM-phenotype', combination='majority')

    ##
    # Export
    ##
    print('Exporting (Complete)')
    wGfile_gpickle = 'results/net-{celltype:s}.gpickle'.format(celltype=celltype)
    ensurePathExists(wGfile_gpickle)
    nx.write_gpickle(G, wGfile_gpickle)
