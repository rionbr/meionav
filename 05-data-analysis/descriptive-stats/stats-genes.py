# coding=utf-8
# Author: Rion B Correia
# Date: April 06, 2020
#
# Description: Information on FPKM genes and how much it maps to protein coding genes
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from tabulate import tabulate


def df2md(df, y_index=False, *args, **kwargs):
    blob = tabulate(df, headers='keys', tablefmt='pipe', *args, **kwargs)
    if not y_index:
        return '\n'.join(['| {}'.format(row.split('|', 2)[-1]) for row in blob.split('\n')])
    return blob


if __name__ == '__main__':

    minLogTPM = 1
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    species = ['HS', 'MM', 'DM']
    ##
    #  CellType
    ##
    r = []
    for celltype in celltypes:

        for specie in species:

            rFPKMfile = "../../02-core_genes/results/FPKM/{specie:s}/{specie:s}-FPKM-{celltype:s}.csv.gz".format(celltype=celltype, specie=specie)
            df = pd.read_csv(rFPKMfile, index_col='id_gene')
            #
            df['log2(TPM)'] = np.log2(df['TPM'] + 1)
            # Remove Duplicates
            df = df.loc[~df.index.duplicated(), :]
            # Only TPM > log(2) & CellType
            df = df.loc[(df['log2(TPM)'] >= minLogTPM), :]

            df.loc[df['biotype'] != 'protein_coding', 'biotype'] = 'other'
            df_pc = df.loc[df['biotype'] == 'protein_coding', :]
            df_si = df.dropna(subset=['id_string'])
            df_si_pc = df_si.loc[df_si['biotype'] == 'protein_coding', :]

            n = df.shape[0]
            n_pc = df_pc.shape[0]
            n_string = df_si.shape[0]
            n_string_pc = df_si_pc.shape[0]

            r.append((specie, celltype, n, n_pc, n_string, n_string_pc))
    print('# Number of genes (LogTPM>=1) with matches to StringDB and that are protein_coding genes\n')

    df_stat = pd.DataFrame(r, columns=['species', 'celltype', '#-genes-log2(TPM)>=1', '#-pc-genes', '#-genes-with-id_string', '#-genes-with-id_string-&-pc'])
    print(df2md(df_stat, floatfmt='.4f'))
    wCSVfile = 'results/stats-genes.csv'
    df_stat.to_csv(wCSVfile)

    #
    # Print Modules
    #
    """
    print('Modules') # modules count start at zero
    #
    modules_infomap = nx.get_node_attributes(G, name='modules-infomap')
    n_modules_infomap = max(modules_infomap.values()) + 1
    # HS
    modules_HS_louvain = nx.get_node_attributes(HSG, name='modules-HS-louvain')
    modules_HS_infomap = nx.get_node_attributes(HSG, name='modules-HS-infomap')
    n_modules_HS_louvain = max(modules_HS_louvain.values()) + 1
    n_modules_HS_infomap = max(modules_HS_infomap.values()) + 1
    # MM
    modules_MM_louvain = nx.get_node_attributes(MMG, name='modules-MM-louvain')
    modules_MM_infomap = nx.get_node_attributes(MMG, name='modules-MM-infomap')
    n_modules_MM_louvain = max(modules_MM_louvain.values()) + 1
    n_modules_MM_infomap = max(modules_MM_infomap.values()) + 1
    # DM
    modules_DM_louvain = nx.get_node_attributes(DMG, name='modules-DM-infomap')
    modules_DM_infomap = nx.get_node_attributes(DMG, name='modules-DM-louvain')
    df_DM_m = pd.DataFrame(data=dict(louvain=modules_DM_louvain, infomap=modules_DM_infomap))
    df_DM_c = pd.DataFrame({
        'nr-mods-louvain': df_DM_m['louvain'].value_counts().value_counts(),
        'nr-mods-infomap': df_DM_m['infomap'].value_counts().value_counts()
        })
    df_DM_c.index.name = 'mod-size'
    n_modules_DM_infomap = max(modules_DM_louvain.values()) + 1
    n_modules_DM_louvain = max(modules_DM_infomap.values()) + 1

    print('# Modules Infomap: {:d}'.format(n_modules_infomap))
    print('# Modules HS Infomap: {:d}'.format(n_modules_HS_infomap))
    print('# Modules HS Louvain: {:d}'.format(n_modules_HS_louvain))
    #
    print('# Modules MM Infomap: {:d}'.format(n_modules_MM_infomap))
    print('# Modules MM Louvain: {:d}'.format(n_modules_MM_louvain))
    #
    print('# Modules DM Infomap: {:d}'.format(n_modules_DM_infomap))
    print('# Modules DM Louvain: {:d}'.format(n_modules_DM_louvain))
    #
    """

    #
    # Correlations
    #

    """
    fertility = nx.get_node_attributes(DMG, 'mean-fert-rate')
    #
    eigen_centrality = nx.eigenvector_centrality(DMG)
    degree_centrality = nx.degree_centrality(DMG)
    #bet_centrality = nx.betweenness_centrality(DMG, k=100)
    page_rank = nx.pagerank(DMG)

    df = pd.DataFrame(data=dict(fertility=fertility, eigen_centrality=eigen_centrality, degree_centrality=degree_centrality, page_rank=page_rank), index=fertility.keys())
    print(df.corr(method='pearson'))
    
    > df.corr(method='pearson')
                       fertility  eigen_centrality  degree_centrality  page_rank
    fertility           1.000000         -0.162134          -0.223829  -0.242392
    eigen_centrality   -0.162134          1.000000           0.876718   0.714450
    degree_centrality  -0.223829          0.876718           1.000000   0.933441
    page_rank          -0.242392          0.714450           0.933441   1.000000
    """

    print('')