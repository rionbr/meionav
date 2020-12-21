# coding=utf-8
# Author: Rion B Correia
# Date: May 19, 2020
#
# Description: Calculates the number of expressed genes (protein-coding) in pairs of species, for all cell types, across different TPM cut-offs.
#
import numpy as np
import pandas as pd
import networkx as nx
from itertools import combinations
from utils import ensurePathExists, get_network_layer


if __name__ == '__main__':

    species = ['HS', 'MM', 'DM']
    celltypes = ['spermatogonia', 'spermatocyte', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    thresholds = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]

    print('Loading Genome Network')
    rGfile_gpickle = '../../04-network/results/network/net-{network:s}.gpickle'.format(network='genome')
    G = nx.read_gpickle(rGfile_gpickle)

    remove_non_protein_coding = [n for n, d in G.nodes(data=True) if d.get('biotype', '') != 'protein_coding']
    G.remove_nodes_from(remove_non_protein_coding)

    """
    print('Separate Layers')
    HSG = get_network_layer(G, 'HS')
    MMG = get_network_layer(G, 'MM')
    DMG = get_network_layer(G, 'DM')

    Gx = {'HS': HSG, 'MM': MMG, 'DM': DMG}
    """

    r = []
    for specie_i, specie_j in combinations(species, 2):

        print("Calculating for species: {specie_i:s} - {specie_j:s}".format(specie_i=specie_i, specie_j=specie_j))

        for celltype in celltypes:

            print("Calculating for celltype: {celltype:s}".format(celltype=celltype))

            rFPKMifile = '../../02-core_genes/results/FPKM/{specie:s}/{specie:s}-FPKM-{celltype:s}.csv.gz'.format(specie=specie_i, celltype=celltype)
            rFPKMjfile = '../../02-core_genes/results/FPKM/{specie:s}/{specie:s}-FPKM-{celltype:s}.csv.gz'.format(specie=specie_j, celltype=celltype)
            #
            df_i = pd.read_csv(rFPKMifile)
            df_j = pd.read_csv(rFPKMjfile)
            #

            for threshold in thresholds:

                genes_i = df_i.loc[((df_i['biotype'] == 'protein_coding') & (df_i['TPM'] >= threshold)), 'id_gene'].tolist()
                genes_j = df_j.loc[((df_j['biotype'] == 'protein_coding') & (df_j['TPM'] >= threshold)), 'id_gene'].tolist()
                #
                genes_ij = genes_i + genes_j

                # Only genes in this modules
                Gtmp = nx.subgraph(G, genes_ij).copy()

                # Remove intra edges
                remove_intra_edges = [(i, j) for i, j, d in Gtmp.edges(data=True) if d.get('type', None) == 'intra']
                Gtmp.remove_edges_from(remove_intra_edges)

                # Remove isolates
                remove_isolates_nodes = list(nx.isolates(Gtmp))
                Gtmp.remove_nodes_from(remove_isolates_nodes)

                a = set(genes_i)
                b = set(genes_j)
                a_union_b = a.union(b)
                a_inter_b = set(Gtmp.nodes())

                similarity = len(a_inter_b) / len(a_union_b)
                r.append((specie_i, specie_j, celltype, threshold, similarity))
    #
    dfR = pd.DataFrame(r, columns=['specie_i', 'specie_j', 'celltype', 'threshold', 'similarity'])

    #
    # Export
    #
    wCSVFile = 'results/pairwise-celltypes-threshold-similarity.csv.gz'
    ensurePathExists(wCSVFile)
    dfR.to_csv(wCSVFile)
