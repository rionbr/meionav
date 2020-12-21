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
#
import argparse


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
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte']  # , 'spermatogonia', 'spermatid']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--species", default='HS', type=str, choices=['HS', 'MM', 'DM'], help="Species to compute GEOA. Defaults to 'HS'.")
    #
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    species = args.species
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

    annotation = '../data/GeneOntology/' + dict_annotation_file[species]
    annotation_reactome = '../data/GeneOntology/reactome.gaf'
    ontology = '../data/GeneOntology/go-basic.obo'

    # for specie in species
    print('Calculating GOEA on {celltype:s} {species:s}'.format(celltype=celltype, species=species))

    # Load Gene Population
    print("Load gene population (from biomart)")
    datamart_name = dict_datamart_names[species]
    ds = Dataset(name=datamart_name, host='http://www.ensembl.org')
    if species == 'DM':
        attributes = ['ensembl_gene_id', 'uniprotswissprot', 'external_gene_name']
    elif species == 'MM':
        attributes = ['ensembl_gene_id', 'uniprotswissprot', 'mgi_id', 'external_gene_name']
    elif species == 'HS':
        attributes = ['ensembl_gene_id', 'uniprotswissprot', 'hmmpanther', 'external_gene_name']
    dfQ = ds.query(attributes=attributes).set_index('Gene stable ID')

    # Population of genes (background) to test against
    if species == 'DM':
        pop_flybase = set(dfQ.index.tolist())
        pop_uniprot = set(dfQ['UniProtKB/Swiss-Prot ID'].dropna().tolist())
        pop = pop_flybase.union(pop_uniprot)
    elif species == 'MM':
        pop_mgi = set(dfQ['MGI ID'].dropna().tolist())
        pop_uniprot = set(dfQ['UniProtKB/Swiss-Prot ID'].dropna().tolist())
        pop = pop_mgi.union(pop_uniprot)
    elif species == 'HS':
        pop_uniprot = set(dfQ['UniProtKB/Swiss-Prot ID'].dropna().tolist())
        pop = pop_uniprot

    # Load GO
    print("Load GO files")
    godag = obo_parser.GODag(ontology)

    # GAF files for both MM and DM
    gaf_species = GafReader(name='GAF ' + species + ' Specie', filename=annotation, godag=godag, namespaces=set(['BP', 'MF']))
    gaf_reactome = GafReader(name='GAF ' + species + ' Reactome', filename=annotation_reactome, godag=godag, namespaces=set(['BP', 'MF']))
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

    """
    # Use PANTHER HMM for HS?
    annotation = '../../data/PANTHER/PANTHER15.0_HMM_classifications'
    dfA = pd.read_csv(annotation, sep='\t', header=None, names=['PANTHER ID', 'name', 'MF', 'BP', 'CC', 'PC', 'PW'], index_col='PANTHER ID')
    dfA = dfA.loc[dfA['BP'].notnull(), ['BP']]

    def keep_only_go(x):
        golist = x.split(';')
        cleanlist = {i.split('#')[-1] for i in golist}
        return cleanlist

    dfA['BP'] = dfA['BP'].apply(keep_only_go)

    dfA['UniProtKB/Swiss-Prot ID'] = dfA.index.map(dfQ.reset_index().set_index('PANTHER ID')['UniProtKB/Swiss-Prot ID'].dropna().to_dict())
    dfA['Gene stable ID'] = dfA.index.map(dfQ.reset_index().set_index('PANTHER ID')['Gene stable ID'].dropna().to_dict())

    ns2assoc_combined = dfA.set_index('UniProtKB/Swiss-Prot ID').loc[:, ['BP']].to_dict()
    """

    # Gene Ontology Enrichment Analysis (GOEA)
    goea = GOEnrichmentStudyNS(pop=pop, ns2assoc=ns2assoc_combined, godag=godag, propagate_counts=True, alpha=0.05, methods=['fdr_bh'])

    # Load PCA
    rCOREFile = 'results/pipeline-{pipeline:s}/{species:s}_meiotic_genes.csv'.format(pipeline='core', species=species)
    #
    wCSVFile = 'results/goea/goea-{celltype:s}-{species:s}-core-genes.csv.gz'.format(celltype=celltype, species=species)

    df_core = pd.read_csv(rCOREFile, index_col=0)

    if species == 'DM':
        genes_flybase = set(df_core.index.tolist())
        genes_uniprot = set(df_core.index.map(dfQ['UniProtKB/Swiss-Prot ID'].dropna().to_dict()).dropna().tolist())
        genes = genes_flybase.union(genes_uniprot)
    elif species == 'MM':
        genes_mgi = set(df_core.index.map(dfQ['MGI ID'].dropna().to_dict()).to_list())
        genes_uniprot = set(df_core.index.map(dfQ['UniProtKB/Swiss-Prot ID'].dropna().to_dict()).dropna().to_list())
        genes = genes_mgi.union(genes_uniprot)
    elif species == 'HS':
        genes = set(df_core.index.map(dfQ['UniProtKB/Swiss-Prot ID'].dropna().to_dict()).dropna().tolist())

    # Run Comparison (only keep GO significant and from 'Biological Process')
    print("> Runnin GOEA test")
    goea_res = goea.run_study(genes, prt=None)

    cols2rem = ['method_flds', 'kws', 'study_items', 'pop_items', 'goterm']
    # transform goea objs for DataFrame format
    res = [{k: v for k, v in i.__dict__.items() if k not in cols2rem} for i in goea_res]
    dfA = pd.DataFrame(res)

    if len(dfA):
        # Index: Biological Process, Significant at 0.01, GO tree depth < 10
        dfS = dfA.loc[(dfA['p_fdr_bh'] <= 0.05), :]
        # Redo Index
        n = len(dfS)
        index = pd.MultiIndex.from_arrays([[celltype] * n, [species] * n], names=('celltype', 'species'))
        dfS.index = index

        print("Exporting")
        ensurePathExists(wCSVFile)
        dfS.to_csv(wCSVFile)
