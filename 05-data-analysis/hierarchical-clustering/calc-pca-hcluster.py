# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads the PCA results (HS, MM & DM) and computes hierarchical clustering
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['figure.titlesize'] = 'medium'
mpl.rcParams['axes.titlesize'] = 'small'
mpl.rcParams['axes.labelsize'] = 'small'
mpl.rcParams['xtick.labelsize'] = 'x-small'
mpl.rcParams['ytick.labelsize'] = 'x-small'
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['hatch.linewidth'] = 0.5
mpl.rcParams['hatch.color'] = '#969696'
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from sklearn.preprocessing import MinMaxScaler
from scipy.interpolate import interp1d
from data import *
import networkx as nx
import scipy.cluster.hierarchy as shc


def add_node(node, is_root=False, min_dist=0.0):
    label = node.id
    node_size = mms.transform([[node.count]])[0][0]
    # Root
    if is_root:
        node_color = 'tab:purple'
    else:
        node_color = specie_node_color
    
    G.add_node(
        node.id,
        label=label,
        node_size=node_size,
        node_color=node_color,
        edgecolors=edgecolors,
        )
    #
    if not node.is_leaf():
        if node.dist > min_dist:
            add_node(node.left, is_root=False, min_dist=min_dist)
            add_node(node.right, is_root=False, min_dist=min_dist)
            G.add_edge(node.id, node.left.id, weight=node.dist)
            G.add_edge(node.id, node.right.id, weight=node.dist)


if __name__ == '__main__':

    celltype = 'spermatocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    layer = 'DM'
    #
    min_dist = 1.5

    colors = {'HS': '#2ca02c', 'MM': '#7f7f7f', 'DM': '#ff7f0e'}

    for celltype in ['spermatocyte', 'enterocyte']:
        for layer in ['HS', 'MM', 'DM']:
            threshold_str = str(threshold).replace('.', 'p')
            print('Computing h-cluster for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
            
            rPCAFile = 'results/pca/{celltype:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            #rHCFile = 'results/hcluster/{celltype:s}/hcluster-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            wIMGFile = 'images/hcluster/{celltype:s}/img-hcluster-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

            specie_node_color = colors[layer]
            df_pca = pd.read_csv(rPCAFile, index_col=0, encoding='utf-8')

            method = 'ward'
            Z = shc.linkage(df_pca[['1c', '2c', '3c', '4c', '5c', '6c', '7c', '8c', '9c']], method)
            dfZ = pd.DataFrame(Z, columns=['idx-i', 'idx-j', 'dist', 'count'], index=pd.RangeIndex(start=0, stop=len(Z), name='iter'))

            dfZ.loc[dfZ['dist'] > 1, :]

            # Scale Count
            mms = MinMaxScaler(feature_range=(3, 125))
            dfZ['count-norm'] = mms.fit_transform(dfZ[['count']])[:, 0]

            root, nodes = shc.to_tree(Z, rd=True)

            G = nx.Graph()
            add_node(root, is_root=True, min_dist=min_dist)
            #
            pos = nx.nx_agraph.graphviz_layout(G, prog='twopi', root=root.id)
            #pos = nx.nx_agraph.graphviz_layout(G, prog='dot', root=root.id)
            #pos = nx.nx_agraph.graphviz_layout(G, prog='neato', root=root.id)
            #pos = nx.nx_agraph.graphviz_layout(G, prog='sfdp', root=root.id)
            #pos = nx.spring_layout(G)
            #pos = nx.circular_layout(G)
            #
            nodes_kwds = {
                'cmap': 'jet',
                'node_size': list(nx.get_node_attributes(G, name='node_size').values()),
                'node_color': list(nx.get_node_attributes(G, name='node_color').values()),
                #'edgecolors': list(nx.get_node_attributes(G, name='edgecolors').values()),
            }
            edges_kwds = {
                'alpha': 1,
                'width': 0.1,
                'edge_color': '#c7c7c7',
                'connectionstyle': 'arc3,rad=0.2'
            }
            #plot_network(G, pos=pos, wIMGFile='images/test-hcluster-network.pdf', **{'nodes': nodes, 'edges': edges})

            print("Plot Network")
            fig, ax = plt.subplots(figsize=(11, 8.5), nrows=1, ncols=1)
            #
            nx.draw_networkx_nodes(G=G, pos=pos, ax=ax, **nodes_kwds)
            nx.draw_networkx_edges(G=G, pos=pos, ax=ax, **edges_kwds)
            s = 'Only leafs with minimum distance of {min_dist:.1f} shown.'.format(min_dist=min_dist)
            ax.annotate(s, xy=(0, 0), xycoords='axes fraction', xytext=(5,5), textcoords='offset points', ha='left', va='bottom')
            #
            plt.tight_layout()
            ensurePathExists(wIMGFile)
            plt.savefig(wIMGFile, dpi=150)
            plt.close()

            """
            print("Plot Dendogram")
            fig, ax = plt.subplots(figsize=(8.5, 9), nrows=1, ncols=1)
            #
            dn = shc.dendrogram(
                Z=Z, ax=ax,
                #truncate_mode='lastp', p=10,
                labels=df_pca['gene'].tolist(), orientation='right', color_threshold=None, leaf_font_size='xx-small')
            #
            plt.tight_layout()
            plt.savefig(wIMGFile, dpi=150)
            plt.close()
            """