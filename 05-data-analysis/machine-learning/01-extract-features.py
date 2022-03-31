# coding=utf-8
# Author: Rion B Correia
# Date: Jan 13, 2020
#
# Description: Reads all available gene information (network, FPKM, DGE, etc) and extracts features for ML.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, ensurePathExists
import argparse
from itertools import product, chain


def ours_or_literature_phenotype(r):
    if pd.notnull(r['Our DM pheno code']):
        return r['Our DM pheno code']
    elif pd.notnull(r['Others DM pheno code']):
        return r['Others DM pheno code']
    else:
        return np.nan


def direct_or_indirect_phenotype(r):
    if pd.notnull(r['direct-phenotype']):
        return r['direct-phenotype']
    elif pd.notnull(r['indirect-phenotype']):
        return 'indirect'
    else:
        return np.nan


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument('--layer', default='DM', type=str, choices=['HS', 'MM', 'DM'], help="Layer/Species.")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    layer = species = args.layer
    layers = ['HS', 'MM', 'DM']
    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    #
    print('Reading {celltype:s}-{network:s}-{threshold:s} Network'.format(celltype=celltype, network=network, threshold=threshold_str))
    path_net = '../../04-network/results/network/{celltype:s}/'.format(celltype=celltype)
    rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Load Multilayer Graph - Extract Layer Graph
    #
    print('Extracting {layer:s} SubGraph'.format(layer=layer))
    Gt = get_network_layer(G, layer)

    #
    # Backbone data 
    #
    print('Reading backbone')
    path_backbone = "../../04-network/results/network-closure/{celltype:s}/".format(celltype=celltype)
    rBfile = path_backbone + "net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle".format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    B = nx.read_gpickle(rBfile)

    is_metric = {(i, j) for i, j, d in B.edges(data=True) if d.get('is_metric') is True}
    Bm = B.edge_subgraph(is_metric).copy()

    is_ultrametric = {(i, j) for i, j, d in B.edges(data=True) if d.get('is_ultrametric') is True}
    Bum = Bm.edge_subgraph(is_ultrametric).copy()

    #
    # (ortho)Backbone data 
    #
    if celltype == 'spermatocyte':
        print('Reading ortho-backbone')
        path_ortho_backbone = "../../04-network/results/network-closure-ortho/{celltype:s}/".format(celltype=celltype)
        rOfile = path_ortho_backbone + "net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle".format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        OB = nx.read_gpickle(rOfile)

        is_metric_ortho = nx.get_edge_attributes(OB, name='is_metric_ortho')
        nx.set_edge_attributes(Gt, name='is_metric_ortho', values=is_metric_ortho)

        is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])
        is_ortho_metric_edges = [(i, j) for i, j, d in OB.edges(data=True) if d.get('is_metric_ortho') == is_metric_ortho_string]
        set_ortho_metric_nodes = set(list(chain(*is_ortho_metric_edges)))
        is_ortho_metric_nodes = {n: n in set_ortho_metric_nodes for n in Gt.nodes()}
        nx.set_node_attributes(Gt, name='is_metric_ortho', values=is_ortho_metric_nodes)

    #
    # Node data to DataFrame
    #
    df = pd.DataFrame.from_dict(dict(Gt.nodes(data=True)), orient='index')

    #
    # Load DGE
    #
    print('Load DEG data')
    path_dge = '../../02-core_genes/results/DE/'
    rfdeg = path_dge + '{species:s}-DE_genes.csv.gz'.format(celltype=celltype, species=species)
    dfdeg = pd.read_csv(rfdeg, index_col=0)
    #
    dfdeg = dfdeg.loc[dfdeg.index.isin(df.index), :]
    # Set DEG variables
    if species == 'DM':
        df['Middle_vs_Apical'] = dfdeg['Middle_vs_Apical']
        df['Middle_vs_Apical'].fillna(False, inplace=True)

        df['Basal_vs_Middle'] = dfdeg['Basal_vs_Middle']
        df['Basal_vs_Middle'].fillna(False, inplace=True)
        #
        df['logFC_MiddleApical'] = dfdeg['logFC_MiddleApical']
        df['logFC_MiddleApical'].fillna(0, inplace=True)
        #
        df['logFC_BasalMiddle'] = dfdeg['logFC_BasalMiddle']
        df['logFC_BasalMiddle'].fillna(0, inplace=True)
    else:
        df['Cyte_vs_Gonia'] = dfdeg['Cyte_vs_Gonia']
        df['Cyte_vs_Gonia'].fillna(False, inplace=True)

        df['Tid_vs_Cyte'] = dfdeg['Tid_vs_Cyte']
        df['Tid_vs_Cyte'].fillna(False, inplace=True)
        #
        df['logFC_CyteGonia'] = dfdeg['logFC_CyteGonia']
        df['logFC_CyteGonia'].fillna(0, inplace=True)
        #
        df['logFC_TidCyte'] = dfdeg['logFC_TidCyte']
        df['logFC_TidCyte'].fillna(0, inplace=True)

    #
    # Load mdlc-mutant DGE
    #
    rMDLCFile = '../../01-diff-gene-exp/results/mdlc/{layer:s}-DGE-mdlc_vs_control.csv'.format(layer=layer)
    dfM = pd.read_csv(rMDLCFile, index_col=0, usecols=['id', 'gene', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'])
    # Filter only DGE significant
    dfMs = dfM.loc[(dfM['logFC'].abs() > 1) & (dfM['FDR'] <= 0.05) & (dfM['logCPM'] >= 1), :].copy()
    dfMs_up = dfMs.loc[(dfMs['logFC'] > 0), :]
    dfMs_dw = dfMs.loc[(dfMs['logFC'] < 0), :]

    def map_up_down(x):
        if x in dfMs_up.index:
            return 'up'
        elif x in dfMs_dw.index:
            return 'down'
        else:
            return 'no-change'
    df['mdlc-mutant-up/down'] = df.index.map(map_up_down)

    df['logFC_mdlc-mutant'] = dfM['logFC']
    df['logFC_mdlc-mutant'].fillna(0, inplace=True)

    #
    # Load mdlc-mutant splicing-defects
    #
    print('Adding mdlc Splicing Defects results')
    rMDLCFile = '../../01-diff-gene-exp/results/mdlc/{layer:s}-IntronRetention-mdlc_vs_control.csv'.format(layer=layer)
    dfI = pd.read_csv(rMDLCFile, index_col=0, usecols=['id', 'gene'])

    df['mdlc-mutant-splidef'] = df.index.map(lambda x: x in dfI.index)

    #
    # Load FPKM
    #
    print('Load FPKM data')
    path_fpkm = '../../02-core_genes/results/FPKM/'
    df_HS_fpkm = pd.read_csv(path_fpkm + 'HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype))
    df_MM_fpkm = pd.read_csv(path_fpkm + 'MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype))
    df_DM_fpkm = pd.read_csv(path_fpkm + 'DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype))

    if species == 'DM':
        dffpkm = df_DM_fpkm.set_index('id_gene')
    elif species == 'MM':
        dffpkm = df_MM_fpkm.set_index('id_gene')
    elif species == 'HS':
        dffpkm = df_HS_fpkm.set_index('id_gene')
    # Only only genes in network.
    #dffpkm = dffpkm.loc[dffpkm.index.isin(df.index), :]

    #
    # Identify conserved genes
    #
    print('Identify Conserved Genes')
    dict_string_gene_HS = df_HS_fpkm.set_index('id_string')['id_gene'].to_dict()
    dict_string_gene_MM = df_MM_fpkm.set_index('id_string')['id_gene'].to_dict()
    dict_string_gene_DM = df_DM_fpkm.set_index('id_string')['id_gene'].to_dict()

    path_meta = '../../02-core_genes/results/meta-genes/'
    dfM = pd.read_csv(path_meta + 'meta-{celltype:s}-genes.csv.gz'.format(celltype=celltype), index_col='id_eggnog', usecols=['id_eggnog', 'id_string_HS', 'id_string_MM', 'id_string_DM'])

    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    dfM['id_gene_HS'] = dfM['id_string_HS'].apply(lambda x: [dict_string_gene_HS[i] for i in x])
    dfM['id_gene_MM'] = dfM['id_string_MM'].apply(lambda x: [dict_string_gene_MM[i] for i in x])
    dfM['id_gene_DM'] = dfM['id_string_DM'].apply(lambda x: [dict_string_gene_DM[i] for i in x])

    dfM = dfM[['id_gene_HS', 'id_gene_MM', 'id_gene_DM']]
    # Only keep meta genes with homologs in at least two species
    #dfM2 = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') >= 2, :]
    # Only keep meta genes with homologs in all three species
    dfM3 = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 3, :]
    #
    dict_conserved = {gene: True for gene in dfM3['id_gene_' + layer].explode().tolist()}
    #
    df['conserved'] = df.index.map(dict_conserved)
    df['conserved'].fillna(False, inplace=True)

    #
    # Load Core genes
    #
    print('Load Core genes')
    path_core = '../../02-core_genes/results/pipeline-core/'
    rfcore = path_core + '{species:s}_meiotic_genes.csv'.format(species=species)
    dfc = pd.read_csv(rfcore, index_col=0)
    #
    core_genes = dfc.index.tolist()
    # Set core
    df.loc[df.index.isin(core_genes), 'core'] = True
    df['core'] = df['core'].fillna(False)

    #
    # Load PCA
    #
    """
    print('Load PCA')
    path_pca = '../../04-network/results/pca/{celltype:s}/{layer:s}/'.format(celltype=celltype, layer=layer)
    rfpcad = path_pca + 'pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    dfpca = pd.read_csv(rfpcad, index_col=0)
    """

    #
    # Load Modules
    #
    print('Load Modules')
    path_modules = '../../05-data-analysis/entropy-based-modules/results/pca-entropy/{celltype:s}/{layer:s}/'.format(celltype=celltype, layer=layer)
    rfmod = path_modules + 'pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-modules.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    dfmod = pd.read_csv(rfmod, index_col=0)

    modcols = []
    for (mid, mname), dfmodt in dfmod.groupby(['module-id', 'module-name']):
        mod_genes = dfmodt.index.tolist()
        modcol = 'module-' + str(mid)
        df[modcol] = df.index.map(lambda x: 1 if x in mod_genes else 0)
        modcols.append(modcol)
    
    #
    # Set Direct Experimental Evidence
    #
    path_evidence = '../../03-screened-data/data/'
    dfe_DM_our = pd.read_csv(path_evidence + 'core_DM_screened_2020-12-23.csv', index_col=0, usecols=['id_gene', 'Our DM pheno code'])
    dfe_DM_our = dfe_DM_our.loc[(dfe_DM_our['Our DM pheno code'].notnull()), 'Our DM pheno code']
    dfe_DM_lit = pd.read_csv(path_evidence + 'literature_DM_2021-01-05.csv', index_col=0, squeeze=True)
    #
    dfe_DM = pd.concat([dfe_DM_our, dfe_DM_lit], axis=1).apply(ours_or_literature_phenotype, axis='columns')
    #
    dfe_MM = pd.read_csv(path_evidence + 'literature_MM_2021-01-05.csv', index_col=0, squeeze=True)
    dfe_HS = pd.read_csv(path_evidence + 'literature_HS_2021-01-05.csv', index_col=0, squeeze=True)
    #
    if species == 'DM':
        df['direct-phenotype'] = dfe_DM
    elif species == 'MM':
        df['direct-phenotype'] = dfe_MM
    elif species == 'HS':
        df['direct-phenotype'] = dfe_HS

    dfMe = pd.DataFrame(index=dfM.index)

    #
    # Set Indirect Experimental Evidence (from homologs)
    #
    dfMe['dir_pheno_HS'] = dfM['id_gene_HS'].map(lambda x: [dfe_HS.get(i) for i in x if i in dfe_HS])
    dfMe['dir_pheno_MM'] = dfM['id_gene_MM'].map(lambda x: [dfe_MM.get(i) for i in x if i in dfe_MM])
    dfMe['dir_pheno_DM'] = dfM['id_gene_DM'].map(lambda x: [dfe_DM.get(i) for i in x if i in dfe_DM])
    #
    # Only keep meta genes with homologs in at least two species
    dfMe1 = dfMe.loc[dfMe.applymap(len).applymap(bool).sum(axis='columns') > 0, :]
    #
    if species == 'DM':
        # add indirect MM and HS
        dfMet = dfMe[['dir_pheno_HS', 'dir_pheno_MM']]
        dfMet = dfMet.loc[dfMet.applymap(len).applymap(bool).sum(axis='columns') > 0, :]
        dfMet['ind-pheno'] = (dfMe['dir_pheno_HS'] + dfMe['dir_pheno_MM']).apply(lambda x: '/'.join(sorted(set(x))))
        dfMet = dfMet['ind-pheno']
        #
        dfMt = dfM.loc[dfMet.index, 'id_gene_DM'].explode().to_frame()
        dfMt['ind-pheno'] = dfMet
        dfMt = dfMt.set_index('id_gene_DM')
        #
        df['indirect-phenotype'] = dfMt['ind-pheno']
    elif species == 'MM':
        # add indirect DM and HS
        dfMet = dfMe[['dir_pheno_HS', 'dir_pheno_DM']]
        dfMet = dfMet.loc[dfMet.applymap(len).applymap(bool).sum(axis='columns') > 0, :]
        dfMet['ind-pheno'] = (dfMe['dir_pheno_HS'] + dfMe['dir_pheno_DM']).apply(lambda x: '/'.join(sorted(set(x))))
        dfMet = dfMet['ind-pheno']
        #
        dfMt = dfM.loc[dfMet.index, 'id_gene_MM'].explode().to_frame()
        dfMt['ind-pheno'] = dfMet
        dfMt = dfMt.set_index('id_gene_MM')
        #
        df['indirect-phenotype'] = dfMt['ind-pheno']
    elif species == 'HS':
        # add indirect DM and MM
        dfMet = dfMe[['dir_pheno_MM', 'dir_pheno_DM']]
        dfMet = dfMet.loc[dfMet.applymap(len).applymap(bool).sum(axis='columns') > 0, :]
        dfMet['ind-pheno'] = (dfMe['dir_pheno_MM'] + dfMe['dir_pheno_DM']).apply(lambda x: '/'.join(sorted(set(x))))
        dfMet = dfMet['ind-pheno']
        #
        dfMt = dfM.loc[dfMet.index, 'id_gene_HS'].explode().to_frame()
        dfMt['ind-pheno'] = dfMet
        dfMt = dfMt.set_index('id_gene_HS')
        #
        df['indirect-phenotype'] = dfMt['ind-pheno']
    #
    df['phenotype'] = df[['direct-phenotype', 'indirect-phenotype']].apply(direct_or_indirect_phenotype, axis='columns')



    #
    # Computing Network Features
    #
    print('Computing Network Features')

    # Degree
    print('> degree')
    dict_degree = {n: d for n, d in list(Gt.degree())}
    df['degree'] = df.index.map(dict_degree)

    print('> degree (metric)')
    dict_degree_metric = {n: d for n, d in list(Bm.degree())}
    df['degree-metric'] = df.index.map(dict_degree_metric)

    # Degree weight
    print('> degree weight')
    dict_weight_degree = {n: d for n, d in list(Gt.degree(weight='weight'))}
    df['degree-weight'] = df.index.map(dict_weight_degree)

    print('> degree weight (metric)')
    dict_weight_degree_metric = {n: d for n, d in list(Bm.degree(weight='weight'))}
    df['degree-weight-metric'] = df.index.map(dict_weight_degree_metric)

    # Degree Centrality
    print('> degree centrality')
    dict_degree_centrality = nx.degree_centrality(Gt)
    df['degree-centrality'] = df.index.map(dict_degree_centrality)

    # Degree Centrality (metric)
    print('> degree centrality (metric)')
    dict_degree_centrality_metric = nx.degree_centrality(Bm)
    df['degree-centrality-metric'] = df.index.map(dict_degree_centrality_metric)


    # Betweeness Centrality
    """
    print('> betweeness centrality (takes a while)')
    dict_betweenness_centrality = nx.betweenness_centrality(Gt, weight='weight')
    df['betweenness-centrality'] = df.index.map(dict_betweenness_centrality)
    """

    # Clustering
    print('> Clustering')
    dict_clustering = nx.clustering(Gt, weight='weight')
    df['clustering'] = df.index.map(dict_clustering)

    print('> Clustering (metric)')
    dict_clustering_metric = nx.clustering(Bm, weight='weight')
    df['clustering-metric'] = df.index.map(dict_clustering_metric)

    # Eigenvector Centrality
    print('> eigenvector centrality')
    dict_eigenvector_centrality = nx.eigenvector_centrality(Gt, weight='weight')
    df['eigenvector-centrality'] = df.index.map(dict_eigenvector_centrality)

    # Eigenvector Centrality
    print('> eigenvector centrality (metric)')
    dict_eigenvector_centrality_metric = nx.eigenvector_centrality(Bm, weight='weight')
    df['eigenvector-centrality-metric'] = df.index.map(dict_eigenvector_centrality_metric)

    # Page rank
    print('> page rank')
    dict_pagerank = nx.pagerank(Gt)
    df['pagerank'] = df.index.map(dict_pagerank)

    # Page rank
    print('> page rank (metric)')

    dict_pagerank_metric = nx.pagerank(Bm)
    df['pagerank-metric'] = df.index.map(dict_pagerank_metric)

    # Average neighbor degree
    print('> avg neighbor degree')
    dict_average_neighbor_degree = nx.average_neighbor_degree(Gt, weight='weight')
    df['avg-neighbor-degree'] = df.index.map(dict_average_neighbor_degree)


    #
    # Edge as column features
    #
    """
    print('Exporting adj.matrix as features')
    dfadj = nx.to_pandas_adjacency(Gt)
    """

    #
    # Edge as column features
    #
    """
    print('Exporting adj.matrix as features')
    dfadjB = nx.to_pandas_adjacency(Bm)
    """

    ##
    # Export
    ##
    print('Saving results to .CSV')
    wMLFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-X-y.csv.gz'.format(celltype=celltype, layer=layer)
    ensurePathExists(wMLFile)
    df.to_csv(wMLFile)

    """
    wAMFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-adj-matrix.csv.gz'.format(celltype=celltype, layer=layer)
    ensurePathExists(wAMFile)
    dfadj.to_csv(wAMFile)

    wAMBFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-adj-matrix-metric.csv.gz'.format(celltype=celltype, layer=layer)
    ensurePathExists(wAMBFile)
    dfadjB.to_csv(wAMBFile)
    """