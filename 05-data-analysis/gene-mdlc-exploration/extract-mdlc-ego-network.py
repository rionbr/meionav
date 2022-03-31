
# coding=utf-8
# Author: Rion B Correia
# Date: June 17, 2021
#
# Description: Extracts mdlc et al subgraph
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

    #
    # Network
    #
    rGfile_gpickle = '../../04-network/results/network/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)
    G = get_network_layer(G, layer)
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


    #
    # mdlc DGE Results
    #
    if layer == 'DM':
        rMDLCFile = '../../01-diff-gene-exp/results/mdlc/{layer:s}-DGE-mdlc_vs_control.csv'.format(layer=layer)
        dfM = pd.read_csv(rMDLCFile, index_col=0, usecols=['id', 'gene', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'])
        # Filter only DGE significant
        dfM = dfM.loc[(dfM['logFC'].abs() > 1) & (dfM['FDR'] <= 0.05) & (dfM['logCPM'] >= 1), :].copy()
        dfM['up/down'] = dfM['logFC'].map(lambda x: 'up' if x > 0 else 'down')
        #
        is_mdlc_dge = dfM['up/down'].to_dict()
        #
        nx.set_node_attributes(G, name='is_mdlc_dge', values=is_mdlc_dge)

    # DEBUG : search for a gene id
    """
    for n,d in Gt.nodes(data=True):
        if d.get('label') == 'Ubi-p63E':
            print(n,d)
            break
    """

    colors = {
        'A': '#2ca02c',  # green as in M2
        'B': '#98df8a',  # lighter green
        #
        'red': '#d62728',
        'other': '#e1e0e2',
    }
    # subgraph
    if layer == 'DM':
        nodes = {
            'FBgn0038772',  # mdlc
            #'FBgn0011206',  # bol
            #'FBgn0003943'  # Ubi-p63E
        }
        node_color = {
            #'mdlc': 'red',
            'dRNF113': 'red',
            # A
            'Usp39': 'A',
            'CG10907': 'A',
            'CG10466': 'A',
            'CG3511': 'A',
            'CG7971': 'A',
            'Cwc25': 'A',
            'CG14641': 'A',
            'c12.1': 'A',
            'CG12343': 'A',
            'Prp18': 'A',
            'CG3225': 'A',
            'CG8435': 'A',
            'pea': 'A',
            'Slu7': 'A',
            'ncm': 'A',
            'Bx42': 'A',
            'Prp8': 'A',
            'fand': 'A',
            'CG4849': 'A',
            'crn': 'A',
            'Dhx15': 'A',
            'l(1)10Bb': 'A',
            'Cdc5': 'A',
            'l(2)37Cb': 'A',
            'CG6015': 'A',
            'l(1)G0007': 'A',
            # B
            'Sf3a1': 'B',
            'Sf3a2': 'B',
            'Sf3b1': 'B',
            'Sf3b2': 'B',
            'Sf3b3': 'B',
            'Sf3b5': 'B',
            'Sf3b6': 'B',
            'Spx': 'B',
            'Phf5a': 'B',
            'SmD1': 'B',
            'SmD2': 'B',
            'SmD3': 'B',
            'SmE': 'B',
            'SmF': 'B',
            'Lsm10': 'B',
            'noi': 'B',
            'l(3)72Ab': 'B',
            'obe': 'B',
            'SNRPG': 'B',
            'scaf6': 'B',
        }
    elif layer == 'MM':
        nodes = {
            'ENSMUSG00000098134'  #
        }
        node_color = {
            'Rnf113a2': 'red',
            # A
            'Cdc5l': 'A',
            'Xab2': 'A',
            'Bud31': 'A',
            'Cdc40': 'A',
            'Cwc15': 'A',
            'Plrg1': 'A',
            'Srrm2': 'A',
            'Rbm22': 'A',
            'Slu7': 'A',
            'Yju2': 'A',
            'Dhx38': 'A',
            'Prpf19': 'A',
            'Eftud3': 'A',
            'Cwc27': 'A',
            'Rbmx2': 'A',
            'Snip1': 'A',
            'Cwc25': 'A',
            'Bud13': 'A',
            'Snw1': 'A',
            'Crnkl1': 'A',
            # B
            'Sf3a2': 'B',
            'Sf3b1': 'B',
            'Sf3b2': 'B',
            'Sf3b3': 'B',
            'Sf3b4': 'B',
            'Sf3b5': 'B',
            'Phf5a': 'B',
            'Snrpf': 'B',
            'Snrpd1': 'B',
            'Snrnp200': 'B',
            'Eftud2': 'B',
        }
    elif layer == 'HS':
        nodes = {
            'ENSG00000139797'
        }
        node_color = {
            #
            'RNF113B': 'red',
            # A
            'CWC22': 'A',
            'SLU7': 'A',
            'YJU2': 'A',
            'DHX38': 'A',
            'BUD13': 'A',
            'SNIP1': 'A',
            'RBMX2': 'A',
            'CWC25': 'A',
            'CWC15': 'A',
            'BUD31': 'A',
            'PRPF19': 'A',
            'CDC40': 'A',
            'SNW1': 'A',
            'CDC5L': 'A',
            'CRNKL1': 'A',
            'PLRG1': 'A',
            'CWC27': 'A',
            'RBM22': 'A',
            # B
            'SF3B1': 'B',
            'SF3B2': 'B',
            'SF3B3': 'B',
            'SF3B4': 'B',
            'PHF5A': 'B',
            'SNRNP200': 'B',
        }


    nodes_and_neighbors = set()
    for ego_node in nodes:
        print(ego_node)
        #
        Gtmp = nx.ego_graph(G, ego_node)
        print(Gtmp.number_of_nodes())
        #
        nodes_and_neighbors.update(Gtmp.nodes())

    #
    # Multilayer Subgraph
    #
    Ge = G.subgraph(nodes_and_neighbors).copy()

    # Manually rename mdlc to dRNF113
    if layer == 'DM':
        Ge.nodes['FBgn0038772']['label'] = 'dRNF113'

    # Remove non backbone edges
    Ge.remove_edges_from([(i, j) for i, j, d in Ge.edges(data=True) if not d.get('is_metric') == True])

    # Remove isolates
    Ge.remove_nodes_from(list(nx.isolates(Ge)))

    #
    # Add Module information to nodes
    #
    """
        Gte = data[layer]
        print('Load/Set Modules ({layer:s})'.format(layer=layer))

        rMfile = '../entropy-based-modules/results/pca-entropy/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-modules.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

        dfm = pd.read_csv(rMfile, index_col=0)
        dfm = dfm.loc[dfm.index.isin(ego_nodes), ['gene', 'module-id']].reset_index()

        dfmg = dfm.groupby('index').agg({'module-id':lambda x: 'module-' + '-'.join(x.astype(str))})

        nx.set_node_attributes(Ge, values=dfmg['module-id'].to_dict(), name='module')
    """
    #
    # Calc is present in other layers
    #
    """
    dict_present = {}
    for node, node_data in Gei.nodes(data=True):
        cur_layer = node_data['layer']
        present = set([cur_layer])
        for neighbor in Gei.neighbors(node):
            other_layer = Gei.nodes[neighbor]['layer']
            present.add(other_layer)
        present_str = 'present-{n:d}'.format(n=len(present))
        dict_present[node] = present_str
    # Set data
    nx.set_node_attributes(Ge, values=dict_present, name='present')
    """

    for i, d in Ge.nodes(data=True):
        is_mdlc_dge = d.get('is_mdlc_dge', None)
        if is_mdlc_dge is not None:
            if is_mdlc_dge == 'up':
                color = '#d62728'
            else:
                color = '#1f77b4'
        else:
            color = colors['other']
        Ge.nodes[i]['color'] = color

    # Manually definedColor nodes
    """
    for i, d in Ge.nodes(data=True): 
        label = d.get('label', None)

        if label in node_color:
            color = colors[node_color[label]]
        else:
            print('{}, {}'.format(i,label))
            color = colors['other']
        Ge.nodes[i]['color'] = color
    """

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
    wGefile_graphml = 'results/net-mdlc-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    ensurePathExists(wGefile_graphml)
    nx.write_graphml(Ge, wGefile_graphml)
