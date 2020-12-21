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
from collections import Counter
#
from nltk.corpus import stopwords
#
from goatools import obo_parser
from goatools.anno.gaf_reader import GafReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS


if __name__ == '__main__':

    celltype = 'spermatocyte'
    network = 'thr'  # 'complete', 'backbone', 'threshold'
    threshold = 0.5
    layer = 'DM'
    threshold_str = str(threshold).replace('.', 'p')
    

    # GO Information
    dict_annotation_file = {'HS': 'goa_human.gaf', 'MM': 'mgi.gaf', 'DM': 'fb.gaf'}
    annotation = '../data/GeneOntology/' + dict_annotation_file[layer]
    ontology = '../data/GeneOntology/go-basic.obo'
    #
    godag = obo_parser.GODag(ontology)
    gaf = GafReader(name='GAF ' + layer, filename=annotation, godag=godag)
    # Dict of Associations
    ns2assoc = gaf.get_ns2assc()

    # Load Population of Genes (for background comparison)
    rFPKMFile = '../02-core_genes/results/FPKM/{layer:s}/{layer:s}-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype, layer=layer)
    dfP = pd.read_csv(rFPKMFile, usecols=['id_gene', 'gene'])

    # Load PCA
    rPCAFile = 'results/pca/{celltype:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    df_pca = pd.read_csv(rPCAFile, index_col=0, encoding='utf-8')

    # Population of genes (background) to test against
    pop = dfP['id_gene'].tolist()

    # Gene Ontology Enrichment Analysis (GOEA)
    goea = GOEnrichmentStudyNS(pop=pop, ns2assoc=ns2assoc, godag=godag, propagate_counts=False, alpha=0.05, methods=['fdr_bh'])

    # Select genes to test
    genes = df_pca.loc[df_pca['1c'] > 5, :].index.tolist()

    # Run Comparison (only keep GO significant and from 'Biological Process')
    goea_res = goea.run_study(genes, prt=None)

    cols2rem = ['method_flds', 'kws', 'study_items', 'pop_items', 'goterm']
    # transform goea objs for DataFrame format
    res = [{k:v for k, v in i.__dict__.items() if k not in cols2rem} for i in goea_res]
    dfA = pd.DataFrame(res)
    # Index: Biological Process, Significant at 0.01, GO tree depth < 10
    dfS = dfA.loc[((dfA['p_fdr_bh'] <= 0.01) & (dfA['depth'] < 10) & (dfA['NS'] == 'BP')), :]
    #
    n = len(dfS)
    index = pd.MultiIndex.from_arrays([[celltype]*n, [layer]*n, [module]*n], names=('celltype', 'layer', 'module'))
    #
    dfS.index = index



    en_stopwords = stopwords.words('english')
    names = [r.name for r in goea_bp_sig]
    wordcount = Counter([word for name in names for word in name.split() if word not in en_stopwords])
    #print('Reading Network')
    #rGfile_gpickle = 'results/net_{network:s}-{level:s}_mlayer.gpickle'.format(network=network, level=levelstr)
    #G = nx.read_gpickle(rGfile_gpickle)


    # Export
    print('Saving results to .CSV')
    #wGOFile = 'results/go/go-{network:s}-{level:s}-{layer:s}.csv.gz'.format(network=network, level=levelstr, layer=layer)
    #ensurePathExists(wGOFile)
    #dfGO.to_csv(wSVDFile)

    print('Done.')