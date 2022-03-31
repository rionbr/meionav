
# coding=utf-8
# Author: Rion B Correia
# Date: March 22, 2021
#
# Description: Reads a MultiLayer module subgraph networks and extracts each layer independently for plotting adding information to color the nodes.
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

    M = G.copy()
    to_remove = {(i, j) for i, j, d in M.edges(data=True) if d.get('type') != 'cross'}
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
        'red': '#d62728',
        #'other': '#e1e0e2', #light gray
        'other': '#7f7f7f', #darker gray
    }
    if layer == 'HS':
        ubi = ['ENSG00000150991', 'ENSG00000170315', 'ENSG00000263563']
        #cdc = ['ENSG00000099804', 'ENSG00000107341']
        #
        node_color = {
            'UBC': 'red',
            'UBB': 'red',
            'UBBP4': 'red',
            'CDC34': 'red',
            'UBE2R2': 'red',
        }
    elif layer == 'MM':
        ubi = ['ENSMUSG00000008348', 'ENSMUSG00000019505']
        #cdc = ['ENSMUSG00000020307', 'ENSMUSG00000020870', 'ENSMUSG00000036241']
        #
        node_color = {
            'Ubc': 'red',
            'Ubb': 'red',
            'Cdc34': 'red',
            'Cdc34b': 'red',
            'Ube2r2': 'red'
        }
    elif layer == 'DM':
        ubi = ['FBgn0003943']
        #cdc = ['FBgn0036516']
        #
        node_color = {
            'Ubi-p63E': 'red',
            #'CG7656': 'red'
        }

    #
    # Ubi
    #
    graphs = []
    for nid in ubi:
        Gts = nx.ego_graph(Gt, nid).copy()  # Ubi
        
        # remove all edges not connected to source (make it a star)
        #to_remove = {(i, j) for i, j in Gts.edges() if not ((i == nid) or (j == nid))}
        #Gts.remove_edges_from(to_remove)

        # remove all edges not metric
        to_remove = {(i, j) for i, j, d in Gts.edges(data=True) if d.get('is_metric', None) != True}
        Gts.remove_edges_from(to_remove)

        # remove all edges not is_metric_ortho
        #is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])
        #to_remove = {(i, j) for i, j, d in Gts.edges(data=True) if d.get('is_metric_ortho', False) != is_metric_ortho_string}
        #Gts.remove_edges_from(to_remove)
        
        # Keep only ortho-neighbors
        #Gts = nx.ego_graph(Gts, nid).copy()
        
        # Giant Connected Component
        conn_comp = sorted(nx.connected_components(Gts), key=len, reverse=True)
        Gts = Gts.subgraph(conn_comp[0]).copy()

        # Remove isolates
        #Gts.remove_nodes_from(list(nx.isolates(Gts)))

        graphs.append(Gts)
    #
    # CG7656
    #
    """
    for nid in cdc:
        Gt_cdc = nx.ego_graph(Gt, cdc).copy()  # CG7656 (CDC34)
        to_remove = {(i, j) for i, j, d in Gt_cdc.edges(data=True) if d.get('is_metric', None) is not True}
        Gt_cdc.remove_edges_from(to_remove)
        # Remove isolates
        Gt_cdc.remove_nodes_from(list(nx.isolates(Gt_cdc)))
    """
    # Compose the two subgraphs
    Gtf = nx.compose_all(graphs)

    # Color nodes
    for i, d in Gtf.nodes(data=True):
        label = d.get('label', None)
        if label in node_color:
            color = colors[node_color[label]]
        else:
            color = colors['other']
        Gtf.nodes[i]['color'] = color

    # Color edges
    is_metric_ortho_string = 'is_metric_ortho' + ''.join(['-{other_layer:s}'.format(other_layer=other_layer) for other_layer in layers if other_layer != layer])
    for i, j, d in Gtf.edges(data=True):
        if d.get('is_metric_ortho', None) == is_metric_ortho_string:
            Gtf[i][j]['color'] = '#e67d7e' # light red  #d62728'  # red
        elif d.get('is_metric', None) == True:
            Gtf[i][j]['color'] = '#80c680' # light green '#2ca02c'  # green
        else:
            Gtf[i][j]['color'] = '#c7c7c7'  # light gray
    
        Gtf[i][j]['weight'] = Gtf[i][j]['weight'] * 10 / 2

    # Export
    wGefile_graphml = 'results/gene-UBI/net-UBI-metric-lcc-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    ensurePathExists(wGefile_graphml)
    nx.write_graphml(Gtf, wGefile_graphml)
