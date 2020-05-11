# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and computed its modules using Louvain & Infomap.
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from utils import ensurePathExists, get_network_layer
#
from goatools import obo_parser
from goatools.anno.gaf_reader import GafReader
from goatools.semantic import TermCounts, get_info_content


if __name__ == '__main__':

    network = 'threshold'  # 'complete', 'backbone', 'threshold'
    level = 0.5
    levelstr = str(level).replace('.', 'p')

    # GO Information
    godag = obo_parser.GODag('../data/GeneOntology/go-basic.obo')
    gaf = GafReader(name='GAF DM', filename='../data/GeneOntology/fb.gaf', godag=godag)

    # dict_gene_go = gaf.get_ns2assc()
    dict_gene_go_nss = gaf.get_id2gos_nss()
    # dict_gene_go_nss['FBgn0037756']
    # term_counts = TermCounts(godag, dict_gene_go_nss)

    print('Reading Network')
    rGfile_gpickle = 'results/net_{network:s}-{level:s}_mlayer.gpickle'.format(network=network, level=levelstr)
    G = nx.read_gpickle(rGfile_gpickle)

    layer = 'DM'
    print('Isolate {layer:s} Layer'.format(layer=layer))
    Gt = get_network_layer(G, layer=layer)
    #
    dfGO = pd.DataFrame(data=None, index=Gt.nodes)

    # Map GO
    dfGO['GO'] = dfG.index.map(lambda x: dict_gene_go_nss.get(x, np.nan))

    # Remove Nan
    dfGO.dropna(subset=['GO'], inplace=True)

    # Export
    print('Saving results to .CSV')
    wGOFile = 'results/go/go-{network:s}-{level:s}-{layer:s}.csv.gz'.format(network=network, level=levelstr, layer=layer)
    ensurePathExists(wGOFile)
    dfGO.to_csv(wSVDFile)

    print('Done.')