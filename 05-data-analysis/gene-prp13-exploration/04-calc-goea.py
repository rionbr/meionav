# coding=utf-8
# Author: Rion B Correia
# Date: Jul 08, 2020
#
# Description: Calculates GO enrichment analysis on PCA modules
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
#
from goatools import obo_parser
from goatools.anno.gaf_reader import GafReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
#
from pybiomart import Dataset


def merge(a, b, path=None):
    # "merges dict b into dict a"
    if path is None:
        path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif isinstance(a[key], set) and isinstance(b[key], set):
                a[key] = a[key].union(b[key])
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a


if __name__ == '__main__':
    #
    # Args
    #
    celltype = 'spermatocyte'
    network = 'rnf113'
    layer = 'HS'
    #
    dict_datamart_names = {
        'HS': 'hsapiens_gene_ensembl',
        'MM': 'mmusculus_gene_ensembl',
        'DM': 'dmelanogaster_gene_ensembl'
    }
    dict_annotation_file = {
        'HS': 'goa_human.gaf',
        'MM': 'mgi.gaf',
        'DM': 'fb.gaf'
    }

    annotation = '../../data/GeneOntology/' + dict_annotation_file[layer]
    annotation_reactome = '../../data/GeneOntology/reactome.gaf'
    ontology = '../../data/GeneOntology/go-basic.obo'

    # for specie in species
    print('Calculating GOEA on orthoBB - DownReg - HS')

    # Load Gene Population
    print("Load gene population (from biomart)")
    datamart_name = dict_datamart_names[layer]
    ds = Dataset(name=datamart_name, host='http://www.ensembl.org')
    attributes = ['ensembl_gene_id', 'uniprotswissprot', 'hmmpanther', 'external_gene_name']
    dfQ = ds.query(attributes=attributes).set_index('Gene stable ID')

    # Population of genes (background) to test against
    pop_uniprot = set(dfQ['UniProtKB/Swiss-Prot ID'].dropna().tolist())
    pop = pop_uniprot

    # Load GO
    print("Load GO files")
    godag = obo_parser.GODag(ontology)

    # GAF files for both MM and DM
    gaf_species = GafReader(name='GAF ' + layer + ' Specie', filename=annotation, godag=godag, namespaces=set(['BP']))
    gaf_reactome = GafReader(name='GAF ' + layer + ' Reactome', filename=annotation_reactome, godag=godag, namespaces=set(['BP']))
    # Dict of Associations
    ns2assoc_species = gaf_species.get_ns2assc()
    n_assoc_species = sum([len(v) for k, v in ns2assoc_species['BP'].items()])
    print('Specie associations: {n:d}'.format(n=n_assoc_species))

    # We also need to add the multi-species annotations
    ns2assoc_reactome = gaf_reactome.get_ns2assc()
    n_assoc_reactome = sum([len(v) for k, v in ns2assoc_reactome['BP'].items()])
    print('Reactome associations: {n:d}'.format(n=n_assoc_reactome))

    # combine associations
    ns2assoc_combined = merge(ns2assoc_species, ns2assoc_reactome)
    n_assoc_combined = sum([len(v) for k, v in ns2assoc_combined['BP'].items()])
    print('Combined associations: {n:d}'.format(n=n_assoc_combined))

    # Gene Ontology Enrichment Analysis (GOEA)
    goea = GOEnrichmentStudyNS(pop=pop, ns2assoc=ns2assoc_combined, godag=godag, propagate_counts=True, alpha=0.05, methods=['fdr_bh'])

    # Load PCA
    rDfile = 'results/net-multilayer-rnf113.csv'
    wCSVFile = 'results/goea-rnf113-HS.csv'

    df = pd.read_csv(rDfile, index_col=0)
    dft = df.loc[(df['layer'] == 'HS') & (df['is_metric_ortho'] == True) & (df['reg'] == 'down') & (df['cross-reg'].str.contains('down')), :]

    genes = set(dft.index.map(dfQ['UniProtKB/Swiss-Prot ID'].dropna().to_dict()).dropna().tolist())

    # Run Comparison (only keep GO significant and from 'Biological Process')
    print("> Runnin GOEA test")
    goea_res = goea.run_study(genes, prt=None)

    cols2rem = ['method_flds', 'kws', 'study_items', 'pop_items', 'goterm']
    # transform goea objs for DataFrame format
    res = [{k: v for k, v in i.__dict__.items() if k not in cols2rem} for i in goea_res]
    dfA = pd.DataFrame(res)

    # Index: Biological Process, Significant at 0.01, GO tree depth < 10
    dfS = dfA.loc[(
        (dfA['p_fdr_bh'] <= 0.05) &
        #(dfA['depth'] < 10) &
        (dfA['NS'] == 'BP')),
    :]
    # Redo Index
    n = len(dfS)

    print("Exporting")
    ensurePathExists(wCSVFile)
    dfS.to_csv(wCSVFile)
