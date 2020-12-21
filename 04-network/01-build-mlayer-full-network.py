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
from utils import ensurePathExists, open_undefined_last_column_files
import argparse


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--minTPM", default=1, type=int, help="minTPM = 1. Defaults to 1.")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    minTPM = args.minTPM
    network = 'full'
    #
    # Init
    #
    G = nx.Graph()

    ##
    # HS Network Data
    ##
    print('Processing HS data')
    df_HS = pd.read_csv("../02-core_genes/results/FPKM/HS/HS-FPKM-{celltype:s}.csv.gz".format(celltype=celltype), index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'FPKM', 'TPM', 'biotype'])
    df_HS = df_HS.rename(columns={'gene': 'label'})  # Rename
    # Only TPM > log(2) & Spermatocyes
    df_HS = df_HS.loc[(df_HS['TPM'] >= minTPM), :]
    # Map id_gene -> id_string
    dict_HS_id_gene_to_id_string = df_HS.explode('id_string').reset_index().set_index('id_string')['id_gene'].to_dict()
    # Bag of String (some genes have >1 id_string)
    set_HS_id_strings = set(np.hstack(df_HS['id_string'].dropna()))
    # logFPKM
    df_HS['logFPKM'] = df_HS['FPKM'].apply(lambda x: np.log2(x + 1))
    df_HS.drop(['FPKM'], axis='columns', inplace=True)
    # Identify layer
    df_HS['layer'] = 'HS'

    # Add Nodes with tuple(id, attrs)
    idxs = df_HS.index.to_list()
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_HS.to_dict(orient='rows')]
    G.add_nodes_from([(i, d) for i, d in zip(idxs, attrs)])

    #
    # HS Links
    #
    print('Adding links')
    df_HS_links = pd.read_csv("../data/StringDB/9606/9606.protein.links.full.v11.0.txt.gz", sep=' ', usecols=['protein1', 'protein2', 'textmining', 'database', 'experiments', 'coexpression', 'neighborhood', 'fusion', 'cooccurence', 'combined_score'])
    # Rename id_string columns
    df_HS_links = df_HS_links.rename(columns={'protein1': 'id_string_i', 'protein2': 'id_string_j'})
    # Reduce Search Space
    df_HS_links = df_HS_links.loc[(df_HS_links['id_string_i'].isin(set_HS_id_strings) & df_HS_links['id_string_j'].isin(set_HS_id_strings)), :]
    # Map id_string to id_gene
    df_HS_links['id_gene_i'] = df_HS_links['id_string_i'].map(lambda x: dict_HS_id_gene_to_id_string[x])
    df_HS_links['id_gene_j'] = df_HS_links['id_string_j'].map(lambda x: dict_HS_id_gene_to_id_string[x])
    #
    df_HS_links.drop(['id_string_i', 'id_string_j'], axis='columns', inplace=True)
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
    df_MM = pd.read_csv("../02-core_genes/results/FPKM/MM/MM-FPKM-{celltype:s}.csv.gz".format(celltype=celltype), index_col='id_gene', usecols=['id_string', 'id_gene', 'gene', 'FPKM', 'TPM', 'biotype'])
    df_MM = df_MM.rename(columns={'gene': 'label'})  # Rename
    # Only TPM >= log(2) & Spermatocyes
    df_MM = df_MM.loc[(df_MM['TPM'] >= minTPM), :]
    # Map id_gene -> id_string
    dict_MM_id_gene_to_id_string = df_MM.explode('id_string').reset_index().set_index('id_string')['id_gene'].to_dict()
    # Bag of String (some genes have >1 id_string)
    set_MM_id_strings = set(np.hstack(df_MM['id_string'].dropna()))
    # logFPKM
    df_MM['logFPKM'] = df_MM['FPKM'].apply(lambda x: np.log2(x + 1))
    df_MM.drop(['FPKM'], axis='columns', inplace=True)
    # Identify layer
    df_MM['layer'] = 'MM'

    # Add Nodes with tuple(id, attrs)
    idxs = df_MM.index.to_list()
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_MM.to_dict(orient='rows')]
    G.add_nodes_from([(i, d) for i, d in zip(idxs, attrs)])

    ##
    # MM Links
    ##
    print('Adding links')
    df_MM_links = pd.read_csv("../data/StringDB/10090/10090.protein.links.full.v11.0.txt.gz", sep=' ', usecols=['protein1', 'protein2', 'textmining', 'database', 'experiments', 'coexpression', 'neighborhood', 'fusion', 'cooccurence', 'combined_score'])
    # Rename id_string columns
    df_MM_links = df_MM_links.rename(columns={'protein1': 'id_string_i', 'protein2': 'id_string_j'})
    # Reduce Search Space
    df_MM_links = df_MM_links.loc[(df_MM_links['id_string_i'].isin(set_MM_id_strings) & df_MM_links['id_string_j'].isin(set_MM_id_strings)), :]
    # Map id_string to id_gene
    df_MM_links['id_gene_i'] = df_MM_links['id_string_i'].map(lambda x: dict_MM_id_gene_to_id_string[x])
    df_MM_links['id_gene_j'] = df_MM_links['id_string_j'].map(lambda x: dict_MM_id_gene_to_id_string[x])
    #
    df_MM_links.drop(['id_string_i', 'id_string_j'], axis='columns', inplace=True)
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
    df_DM = pd.read_csv("../02-core_genes/results/FPKM/DM/DM-FPKM-{celltype:s}.csv.gz".format(celltype=celltype), index_col='id_gene', usecols=['id_gene', 'id_string', 'gene', 'FPKM', 'TPM', 'biotype'])
    df_DM = df_DM.rename(columns={'gene': 'label'})  # Rename
    # Only TPM >= log(2)
    df_DM = df_DM.loc[(df_DM['TPM'] >= minTPM), :]
    # Map id_gene -> id_string
    dict_DM_id_gene_to_id_string = df_DM.explode('id_string').reset_index().set_index('id_string')['id_gene'].to_dict()
    # Bag of String (some genes have >1 id_string)
    set_DM_id_strings = set(np.hstack(df_DM['id_string'].dropna()))
    # logFPKM
    df_DM['logFPKM'] = df_DM['FPKM'].apply(lambda x: np.log2(x + 1))
    df_DM.drop(['FPKM'], axis='columns', inplace=True)
    # Identify layer
    df_DM['layer'] = 'DM'

    # Add Nodes with tuple(id, attrs)
    idxs = df_DM.index.to_list()
    attrs = [{k: v for k, v in m.items() if pd.notnull(v)} for m in df_DM.to_dict(orient='rows')]
    G.add_nodes_from([(i, d) for i, d in zip(idxs, attrs)])

    ##
    # DS Links
    ##
    print('Adding links')
    df_DM_links = pd.read_csv("../data/StringDB/7227/7227.protein.links.full.v11.0.txt.gz", sep=' ', usecols=['protein1', 'protein2', 'textmining', 'database', 'experiments', 'coexpression', 'neighborhood', 'fusion', 'cooccurence', 'combined_score'])
    # Rename id_string columns
    df_DM_links = df_DM_links.rename(columns={'protein1': 'id_string_i', 'protein2': 'id_string_j'})
    # Reduce Search Space
    df_DM_links = df_DM_links.loc[(df_DM_links['id_string_i'].isin(set_DM_id_strings) & df_DM_links['id_string_j'].isin(set_DM_id_strings)), :]
    # Map id_string to id_gene
    df_DM_links['id_gene_i'] = df_DM_links['id_string_i'].map(lambda x: dict_DM_id_gene_to_id_string[x])
    df_DM_links['id_gene_j'] = df_DM_links['id_string_j'].map(lambda x: dict_DM_id_gene_to_id_string[x])
    #
    df_DM_links.drop(['id_string_i', 'id_string_j'], axis='columns', inplace=True)
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
    # dfM = pd.read_csv('../02-core_genes/results/meta_meiotic_genes.csv.gz', index_col='id_eggnog')
    set_id_strings = set.union(set_HS_id_strings, set_MM_id_strings, set_DM_id_strings)

    # Load EggNOG Annotation File
    # Metazoa (33208) EggNOG - [M]embers
    df_Egg = open_undefined_last_column_files(
        "../data/EggNOG/33208_members.tsv.gz",
        n_fixed_cols=5,
        names=['family', 'id_eggnog', '_1', '_2', 'aliases', 'species'],
        nrows=None
    )
    # Only keep columns we need
    df_Egg = df_Egg.set_index('id_eggnog')['aliases']
    # Only keep species we need")
    wanted_species = frozenset(['7227', '9606', '10090'])

    def select_by_species(text, keeplist):
        # Only keep genes from species we are interested in (lower the search space)
        return [i for i in text.split(',') if i.split('.', 1)[0] in keeplist]

    df_Egg = df_Egg.apply(select_by_species, args=(wanted_species,))

    #  Separating by At Least One Match
    def select_by_at_least_one_match(ilist, keeplist):
        # Only keep genes that are found in any of our gene list (lower the search space)
        spermgenes = [i for i in ilist if i in keeplist]
        return spermgenes if len(spermgenes) >= 1 else None

    df_Egg = df_Egg.apply(select_by_at_least_one_match, args=(set_id_strings, ))
    df_Egg = df_Egg.dropna()

    # Selecting by species
    def select_by_gene_and_separate_by_species(ilist, keeplist_HS, keeplist_MM, keeplist_DM):
        # Separate by species, keeping only the sperm genes we are interested in
        genes_HS = [i for i in ilist if i in keeplist_HS]
        genes_MM = [i for i in ilist if i in keeplist_MM]
        genes_DM = [i for i in ilist if i in keeplist_DM]

        return pd.Series({'id_string_HS': genes_HS, 'id_string_MM': genes_MM, 'id_string_DM': genes_DM})

    dfM = df_Egg.apply(select_by_gene_and_separate_by_species, args=(set_HS_id_strings, set_MM_id_strings, set_DM_id_strings))

    # Add the edges
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
    # Export
    ##
    print('Exporting')
    wGfile_gpickle = 'results/network/{celltype:s}/net-{celltype:s}-{network:s}.gpickle'.format(celltype=celltype, network=network)
    ensurePathExists(wGfile_gpickle)
    nx.write_gpickle(G, wGfile_gpickle)

    print('Done.')
