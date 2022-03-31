# coding=utf-8
# Author: Rion B Correia
# Date: Feb 03, 2022
#
# Description: Plots an overview of the FPKM across the three species
#
#
import numpy as np
import pandas as pd
pd.set_option("display.max_columns", 10)
pd.set_option("display.width", 150)
import networkx as nx
from itertools import product
#import matplotlib as mpl
#mpl.rcParams['font.family'] = 'Helvetica'
#mpl.rcParams['mathtext.fontset'] = 'cm'
#import matplotlib.pyplot as plt


if __name__ == '__main__':

    layers = ['HS', 'MM', 'DM']
    layer = 'DM'

    df = pd.read_csv('data/selected-ortho-edges.csv', index_col=0)

    dfMHS = pd.read_csv('../gene-prp13-exploration/results/HS-DGE-rnf113_vs_control-SELECTED-string.csv', index_col=0)

    df['gene_i-rnf113-dge'] = df['id_gene_i'].isin(dfMHS.index)
    df['gene_j-rnf113-dge'] = df['id_gene_j'].isin(dfMHS.index)

    #
    # List of HS that appeared in RNF113-DGE
    #
    df.loc[(df['gene_i-rnf113-dge'] == True) | (df['gene_j-rnf113-dge'] == True), :].to_csv('results/selected-ortho-edges-HS-rnf113.csv')
    #df.loc[(df['gene_i-rnf113-dge'] == True) | (df['gene_j-rnf113-dge'] == True), ['id_gene_i', 'id_gene_j', 'gene_i', 'gene_j', 'gene_i-rnf113-dge', 'gene_j-rnf113-dge']]


    #
    # Now find the respective selected-ortho-edges in DM
    #
    #
    # Threshold Network
    #
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype='spermatocyte', network='thr', threshold='0p5')
    G = nx.read_gpickle(rGfile_gpickle)

    # OrthoBb
    rOfile_gpickle = '../../04-network/results/network-closure-ortho/{celltype:s}/net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype='spermatocyte', network='thr', threshold='0p5', layer='DM')
    GO = nx.read_gpickle(rOfile_gpickle)

    is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])
    edges2remove = [(i, j) for i, j, d in GO.edges(data=True) if d.get('is_metric_ortho') != is_metric_ortho_string]
    GO.remove_edges_from(edges2remove)

    # Intra layers
    is_not_cross_edges = [(i, j) for i, j, d in G.edges(data=True) if d.get('type') != 'cross']
    H = G.copy()
    H.remove_edges_from(is_not_cross_edges)

    # Loop every edge
    selected_DM_ortho_edges = []
    for idx, r in df.iterrows():
        hs_i = r['id_gene_i']
        hs_j = r['id_gene_j']
        hs_i_name = r['gene_i']
        hs_j_name = r['gene_j']

        dm_is = [n for n in H[hs_i] if H.nodes[n].get('layer') == 'DM']
        dm_js = [n for n in H[hs_j] if H.nodes[n].get('layer') == 'DM']

        for dm_i, dm_j in list(product(dm_is, dm_js)):
            if GO.has_edge(dm_i, dm_j):

                iname = GO.nodes[dm_i].get('label')
                jname = GO.nodes[dm_j].get('label')
                edge = GO[dm_i][dm_j]
                weight = edge['weight']

                selected_DM_ortho_edges.append((dm_i, dm_j, iname, jname, weight, hs_i, hs_j, hs_i_name, hs_j_name))

    dfS = pd.DataFrame(selected_DM_ortho_edges, columns=['id_gene_i', 'id_gene_j', 'gene_i', 'gene_j', 'weight', 'HS_ortho-id_gene_i', 'HS_ortho-id_gene_i', 'HS_ortho-gene_i', 'HS_ortho-gene_j'])
    dfS.to_csv('data/selected-ortho-genes-DM.csv')
    #
    # List of DM taht appear in dRNF113-DGE
    #
    dfMDM = pd.read_csv('../gene-prp13-exploration/results/DM-DGE-mdlc_vs_control-SELECTED-string.csv', index_col=0)

    dfS['gene_i-mdlc-dge'] = dfS['id_gene_i'].isin(dfMDM.index)
    dfS['gene_j-mdlc-dge'] = dfS['id_gene_j'].isin(dfMDM.index)

    dfS = dfS.loc[(dfS['gene_i-mdlc-dge'] == True) | (dfS['gene_j-mdlc-dge'] == True), :]

    dfS.to_csv('results/selected-ortho-genes-DM-mdlc.csv')