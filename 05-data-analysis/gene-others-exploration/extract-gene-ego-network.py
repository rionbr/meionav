
# coding=utf-8
# Author: Rion B Correia
# Date: March 22, 2021
#
# Description: Extracts an ego-network from any of the species graph
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer, ensurePathExists


if __name__ == '__main__':

    celltype = 'spermatocyte'
    network = 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    layer = 'DM'
    layers = ['HS', 'MM', 'DM']
    whichgenes = 'UBI'  #  'CDYL', 'SPATA5' 'OSBP2'
    #
    # Network
    #
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    #
    # Add Backbone
    #
    path_backbone = '../../04-network/results/network-closure/{celltype:}/'.format(celltype=celltype)
    rBfile_gpickle = path_backbone + 'net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    B = nx.read_gpickle(rBfile_gpickle)
    #
    is_metric = nx.get_edge_attributes(B, 'is_metric')
    is_ultrametric = nx.get_edge_attributes(B, 'is_ultrametric')
    nx.set_edge_attributes(G, values=is_metric, name='is_metric')
    nx.set_edge_attributes(G, values=is_ultrametric, name='is_ultrametric')

    #
    # Add Ortho-backbone
    #
    path_ortho_backbone = '../../04-network/results/network-closure-ortho/{celltype:s}/'.format(celltype=celltype)
    rOfile_gpickle = path_ortho_backbone + 'net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    O = nx.read_gpickle(rOfile_gpickle)
    #
    is_metric_ortho = nx.get_edge_attributes(O, 'is_metric_ortho')
    nx.set_edge_attributes(G, values=is_metric_ortho, name='is_metric_ortho')

    # only HS layer
    Gt = get_network_layer(G, layer)


    # DEBUG : search for a gene id
    """
    for n,d in G.nodes(data=True):
        if d.get('label') == 'OSBP2':
            print(n,d)
            break
    """

    # General Colors
    colors = {
        # CDYL
        'A': '#36D7CC',  # Transcriptional repressors
        'B': '#E04CE3',  # Histone proteins
        # SPATA5
        'C': '#FCD700',  # Proteasome subunits
        'D': '#EB7600',  # Ubiquitin regulation
        'E': '#E95F5E',  # Protein degradation
        # OSBP2
        'F': '#f7b6d2',  # Steroid degradation
        'G': '#dbdb8d',  # Protein folding
        'H': '#b9b5ff',  # Membrane trafficking
        #
        'red': '#d62728',
        #
        'light-green': '#eaf4e4',
        'light-red': '#ffd8d4',
        'other': '#e1e0e2',
    }

    if whichgenes == 'CDYL':
        nodes = {
            'ENSG00000153046',  # 'CDYL'
        }
        node_color = {
            'CDYL': 'red',
            # Transcriptional repressors
            'MIER1': 'A',
            'EZH2': 'A',
            'SUZ12': 'A',
            #
            'EHMT1': 'A',
            'SETDB1': 'A',
            'HDAC1': 'A',
            'RCOR1': 'A',
            'REST': 'A',
            # Histone proteins
            'HIST1H4C': 'B',
            'HIST1H4B': 'B',
            'HIST1H4E': 'B',
            'HIST1H4D': 'B',
            'HIST1H4I': 'B',
            'HIST1H3A': 'B',
            'HIST2H4B': 'B',
            'HIST1H4H': 'B',
            'HIST1H3D': 'B',
            #
            'WIZ': 'light-green',
            'ZNF280B': 'light-green',
            'RPL23A': 'light-green',
            'LARP1B': 'light-green',
            'UBAP2': 'light-green',
            'SMARCA4': 'light-green',
            'DNMT1': 'light-green',
        }

    elif whichgenes == 'SPATA5':
        nodes = {
            'ENSG00000145375',  # SPATA5
        }
        node_color = {
            'SPATA5': 'red',
            # Proteasome subunits
            'PSMD1': 'C',
            'PSMD2': 'C',
            'PSMC1': 'C',
            'PSMC2': 'C',
            'PSMC3': 'C',
            'PSMC6': 'C',
            # Ubiquitin regulation
            'TRIP12': 'D',
            'UBE4A': 'D',
            'UBOX5': 'D',
            'UFD1': 'D',
            'NPLOC4': 'D',
            #
            'AFG1L': 'E',
            #(REMOVED)'UBXN6': 'E',
            #
            'TMEM33': 'light-red',
            'TAPT1': 'light-red',
            #
            'CDC5L': 'light-green',
            'CINP': 'light-green',
            'ANKRD50': 'light-green',
            'SPATA5L1': 'light-green',
            'ARL9': 'light-green',
            'LIMCH1': 'light-green',
            'UBXN11': 'light-green',
        }
    elif whichgenes == 'OSBP2':
        nodes = {
            'ENSG00000184792',  # OSBP2
        }
        node_color = {
            'OSBP2': 'red',
            # Steroid degradation
            'UGT1A6': 'F',
            'CYP2C8': 'F',
            # Protein folding
            'PPIG': 'G',
            #(REMOVED)'SELENOF': 'G',
            # Membrane trafficking
            'VAPA': 'H'
        }
    elif whichgenes == 'UBI':
        nodes = {
            'FBgn0003943',  # Ubi-p63E
        }
        node_color = {
            'Ubi-p63E': 'red',
        }

    """
    elif whichgenes == 'DAZAP1':
        nodes = {
            'ENSG00000071626',  # 'DAZAP1'
        }

    elif whichgenes == 'BOLL':
        nodes = {
            'ENSG00000152430',  # 'BOLL'
        }
    elif whichgenes == 'RSBN1+SPATA20':
        # RSBN1 + SPATA20
        nodes = {
            'ENSG00000081019',  # RSBN1
            'ENSG00000006282',  # SPATA20
        }
    elif whichgenes == 'TEX2':
        nodes = {
            'ENSG00000136478',  # TEX2
            #'ENSG00000104450',  # SPAG1
        }
    elif whichgenes == 'FAM50A':
        nodes = {
            'ENSG00000071859',  # FAM50A
        }
    elif whichgenes == 'MYCBPAP':
        nodes = {
            'ENSG00000136449',  # MYCBPAP
        }
    """

    nodes_and_neighbors = set()
    for ego_node in nodes:
        #
        Gtmp = nx.ego_graph(Gt, ego_node)
        #
        nodes_and_neighbors.update(Gtmp.nodes())

    #
    # Multilayer Subgraph
    #
    Ge = Gt.subgraph(nodes_and_neighbors)

    #
    # Add conserved information
    #
    """
    path_fpkm = '../../02-core_genes/results/FPKM/'
    df_HS = pd.read_csv(path_fpkm + 'HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
    df_MM = pd.read_csv(path_fpkm + 'MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
    df_DM = pd.read_csv(path_fpkm + 'DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')

    dict_string_gene_HS = df_HS['id_gene'].to_dict()
    dict_string_gene_MM = df_MM['id_gene'].to_dict()
    dict_string_gene_DM = df_DM['id_gene'].to_dict()

    print('Loading {celltype:s} meta genes'.format(celltype=celltype))
    path = '../../02-core_genes/results/'
    dfM = pd.read_csv(path + 'meta-genes/meta-{celltype:s}-genes.csv.gz'.format(celltype=celltype), index_col='id_eggnog', usecols=['id_eggnog', 'id_string_HS', 'id_string_MM', 'id_string_DM'])

    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    dfM['id_gene_HS'] = dfM['id_string_HS'].apply(lambda x: [dict_string_gene_HS[i] for i in x])
    dfM['id_gene_MM'] = dfM['id_string_MM'].apply(lambda x: [dict_string_gene_MM[i] for i in x])
    dfM['id_gene_DM'] = dfM['id_string_DM'].apply(lambda x: [dict_string_gene_DM[i] for i in x])

    dfM = dfM[['id_gene_HS', 'id_gene_MM', 'id_gene_DM']]
    # Only keep meta genes with homologs in all three species
    dfM = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 3]

    dict_conserved = {gene: True for gene in dfM['id_gene_' + layer].explode().tolist()}
    nx.set_node_attributes(Ge, values=dict_conserved, name='conserved')
    """
    # Color nodes
    for i, d in Ge.nodes(data=True):
        label = d.get('label', None)
        if label in node_color:
            color = colors[node_color[label]]
        else:
            color = colors['other']
        Ge.nodes[i]['color'] = color

    # Color edges
    is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])
    for i, j, d in Ge.edges(data=True):
        if d.get('is_metric_ortho', None) == is_metric_ortho_string:
            Ge[i][j]['color'] = '#d62728'  # red
        elif d.get('is_metric', None) == True:
            Ge[i][j]['color'] = '#2ca02c'  # green
        else:
            Ge[i][j]['color'] = '#c7c7c7'  # light gray

    # Export
    wGefile_graphml = 'results/gene-{whichgenes:s}/net-{whichgenes:s}-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(whichgenes=whichgenes, celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    ensurePathExists(wGefile_graphml)
    nx.write_graphml(Ge, wGefile_graphml)
