# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and prints information.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import get_network_layer
from tabulate import tabulate


def df2md(df, y_index=False, *args, **kwargs):
    blob = tabulate(df, headers='keys', tablefmt='pipe', *args, **kwargs)
    if not y_index:
        return '\n'.join(['| {}'.format(row.split('|', 2)[-1]) for row in blob.split('\n')])
    return blob


if __name__ == '__main__':
    
    celltypes = ['enterocyte', 'spermatocyte']
    #
    # Complete Network
    #
    r = []
    for celltype in celltypes:

        print('Loading {celltype:s} network'.format(celltype=celltype))
        rGfile_gpickle = 'results/net-{celltype:s}.gpickle'.format(celltype=celltype)
        G = nx.read_gpickle(rGfile_gpickle)

        for layer in ['HS', 'MM', 'DM']:
            print('Separate {layer:s} layer'.format(layer=layer))
            Gt = get_network_layer(G, layer)

            # Number of nodes/edges
            n_nodes = Gt.number_of_nodes()
            n_edges = Gt.number_of_edges()

            r.append((celltype, layer, n_nodes, n_edges))

    print('# Number of nodes/edges in each layer of the complete network\n')
    df_stat = pd.DataFrame(r, columns=['Celltype', 'Species', '# Nodes', '# Edges'])
    print(df2md(df_stat, floatfmt='.4f'))

    #
    # Threshold Network
    #
    r = []
    for celltype in celltypes:

        network = 'thr'  # 'thr'
        threshold = 0.5
        threshold_str = str(threshold).replace('.', 'p')

        print('Reading {celltype:s}-{network:s}-{threshold:s} network'.format(celltype=celltype, network=network, threshold=threshold_str))
        rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        G = nx.read_gpickle(rGfile_gpickle)

        for layer in ['HS', 'MM', 'DM']:
            print('Separate {layer:s} layer'.format(layer=layer))
            Gt = get_network_layer(G, layer)

            # Number of nodes/edges
            n_nodes = Gt.number_of_nodes()
            n_edges = Gt.number_of_edges()
            # Number of components (islands)
            n_components = nx.number_connected_components(Gt)
            # Largest Component
            Gtlc = max(nx.connected_component_subgraphs(Gt), key=len)
            n_nodes_largest_component = Gtlc.number_of_nodes()
            n_edges_largest_component = Gtlc.number_of_edges()

            r.append((celltype, layer, n_nodes, n_edges, n_components, n_nodes_largest_component, n_edges_largest_component))

    print('# Number of nodes/edges in the layer of the thresholded>0.5 network\n')

    df_stat = pd.DataFrame(r, columns=['Celltype', 'Species','# Nodes','# Edges', '# Components', '# Nodes Largest Comp.', '#Edges Largest Comp.'])
    print(df2md(df_stat, floatfmt='.4f'))

    #
    # Number of Experimental Nodes/Edges
    #
    r = []
    for celltype in celltypes:

        network = 'thr'  # 'thr'
        threshold = 0.5
        threshold_str = str(threshold).replace('.', 'p')

        print('Reading {celltype:s}-{network:s}-{threshold:s} network'.format(celltype=celltype, network=network, threshold=threshold_str))
        rGfile_gpickle = 'results/net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
        G = nx.read_gpickle(rGfile_gpickle)

        for layer in ['HS', 'MM', 'DM']:
            print('Separate {layer:s} layer'.format(layer=layer))
            Gt = get_network_layer(G, layer)
            Gtexp = Gt.copy()
            Gtexp0p5 = Gt.copy()
            
            edges_to_remove = [(i, j) for i, j, d in Gtexp.edges(data=True) if 'experiments' not in d]
            Gtexp.remove_edges_from(edges_to_remove)
            Gtexp.remove_nodes_from(list(nx.isolates(Gtexp)))

            edges_to_remove = [(i, j) for i, j, d in Gtexp0p5.edges(data=True) if d.get('experiments', 0) <= 500]
            Gtexp0p5.remove_edges_from(edges_to_remove)
            Gtexp0p5.remove_nodes_from(list(nx.isolates(Gtexp0p5)))

            # Number of nodes/edges
            n_nodes_exp = Gtexp.number_of_nodes()
            n_edges_exp = Gtexp.number_of_edges()

            # Number of nodes/edges > 0.5
            n_nodes_exp_0p5 = Gtexp0p5.number_of_nodes()
            n_edges_exp_0p5 = Gtexp0p5.number_of_edges()
            # Number of components (islands)
            #n_components = nx.number_connected_components(Gt)
            # Largest Component
            #Gtlc = max(nx.connected_component_subgraphs(Gt), key=len)
            #n_nodes_largest_component = Gtlc.number_of_nodes()
            #n_edges_largest_component = Gtlc.number_of_edges()

            r.append((celltype, layer, n_nodes_exp, n_edges_exp, n_nodes_exp_0p5, n_edges_exp_0p5))

    print('# Number of experimental nodes/edges of the thresholded>0.5 network\n')
    df_stat = pd.DataFrame(r, columns=['Celltype', 'Species', '# Exp. Nodes', '# Exp. Edges', '# Exp. Nodes (>0.5)', '# Exp. Edges (>0.5)'])
    print(df2md(df_stat, floatfmt='.4f'))
 
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