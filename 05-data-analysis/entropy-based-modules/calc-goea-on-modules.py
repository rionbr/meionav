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
from data_spermatocyte_pca_modules_dm import spermatocyte_pca_modules_dm
from data_spermatocyte_pca_modules_mm import spermatocyte_pca_modules_mm
from data_spermatocyte_pca_modules_hs import spermatocyte_pca_modules_hs
#
from data_enterocyte_pca_modules_dm import enterocyte_pca_modules_dm
from data_enterocyte_pca_modules_mm import enterocyte_pca_modules_mm
from data_enterocyte_pca_modules_hs import enterocyte_pca_modules_hs
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
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    parser.add_argument("--layer", default='HS', type=str, choices=['DM', 'MM', 'HS'], help="Network layer to compute SVD. Defaults to 'DM'.")
    #
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    layer = args.layer
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

    data = {
        'spermatocyte': {
            'HS': spermatocyte_pca_modules_hs,
            'MM': spermatocyte_pca_modules_mm,
            'DM': spermatocyte_pca_modules_dm,
        },
        'enterocyte': {
            'HS': enterocyte_pca_modules_hs,
            'MM': enterocyte_pca_modules_mm,
            'DM': enterocyte_pca_modules_dm,
        }
    }
    modules = data[celltype][layer]

    # for specie in species
    print('Calculating GOEA on {celltype:s} {network:s} {threshold:.1f} {layer:s}'.format(celltype=celltype, network=network, threshold=threshold, layer=layer))

    # Load Gene Population
    print("Load gene population (from biomart)")
    datamart_name = dict_datamart_names[layer]
    ds = Dataset(name=datamart_name, host='http://www.ensembl.org')
    if layer == 'DM':
        attributes = ['ensembl_gene_id', 'uniprotswissprot', 'external_gene_name']
    elif layer == 'MM':
        attributes = ['ensembl_gene_id', 'uniprotswissprot', 'mgi_id', 'external_gene_name']
    elif layer == 'HS':
        attributes = ['ensembl_gene_id', 'uniprotswissprot', 'hmmpanther', 'external_gene_name']
    dfQ = ds.query(attributes=attributes).set_index('Gene stable ID')

    # Population of genes (background) to test against
    if layer == 'DM':
        pop_flybase = set(dfQ.index.tolist())
        pop_uniprot = set(dfQ['UniProtKB/Swiss-Prot ID'].dropna().tolist())
        pop = pop_flybase.union(pop_uniprot)
    elif layer == 'MM':
        pop_mgi = set(dfQ['MGI ID'].dropna().tolist())
        pop_uniprot = set(dfQ['UniProtKB/Swiss-Prot ID'].dropna().tolist())
        pop = pop_mgi.union(pop_uniprot)
    elif layer == 'HS':
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
    rPCAFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rDiAnFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dian.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rEntFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-entropy.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    #
    wCSVFile = 'results/goea/{celltype:s}/goea-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

    df_pca = pd.read_csv(rPCAFile, index_col=0)
    df_dian = pd.read_csv(rDiAnFile, index_col=0)
    df_ent = pd.read_csv(rEntFile, index_col=0)

    ldfS = []
    for module in modules:

        mid = module['id']
        mname = module['name']

        print("Computing Module {mid:d}".format(mid=mid))
        #
        xc = module['dim-coords']['xdim']
        yc = module['dim-coords']['ydim']
        ld1 = '{xc:d}c'.format(xc=xc)  # label 1 component
        ld2 = '{yc:d}c'.format(yc=yc)  # label 2 component
        l12d = '{xc:d}c-{yc:d}c-dist'.format(xc=xc, yc=yc)  # label-1c-2c-dist

        x0, x1 = module['dim-coords']['xvals']
        y0, y1 = module['dim-coords']['yvals']
        cut_rank = module['dim-coords']['radius-rank']

        # Radius of the circle
        cut_radius = df_ent.loc[((df_ent['dim'] == xc) & (df_ent['cut-rank'] == cut_rank)), 'radius-start'].squeeze()

        # Select genes in module
        df_pca_tmp = df_pca.loc[((df_pca[ld1] > x0) & (df_pca[ld1] < x1) & (df_pca[ld2] > y0) & (df_pca[ld2] < y1) & (df_dian[l12d] > cut_radius)), :]
        if layer == 'DM':
            genes_flybase = set(df_pca_tmp.index.tolist())
            genes_uniprot = set(df_pca_tmp.index.map(dfQ['UniProtKB/Swiss-Prot ID'].dropna().to_dict()).dropna().tolist())
            genes = genes_flybase.union(genes_uniprot)
        elif layer == 'MM':
            genes_mgi = set(df_pca_tmp.index.map(dfQ['MGI ID'].dropna().to_dict()).to_list())
            genes_uniprot = set(df_pca_tmp.index.map(dfQ['UniProtKB/Swiss-Prot ID'].dropna().to_dict()).dropna().to_list())
            genes = genes_mgi.union(genes_uniprot)
        elif layer == 'HS':
            genes = set(df_pca_tmp.index.map(dfQ['UniProtKB/Swiss-Prot ID'].dropna().to_dict()).dropna().tolist())

        # Run Comparison (only keep GO significant and from 'Biological Process')
        print("> Runnin GOEA test")
        goea_res = goea.run_study(genes, prt=None)

        cols2rem = ['method_flds', 'kws', 'study_items', 'pop_items', 'goterm']
        # transform goea objs for DataFrame format
        res = [{k: v for k, v in i.__dict__.items() if k not in cols2rem} for i in goea_res]
        dfA = pd.DataFrame(res)

        if len(dfA):
            # Index: Biological Process, Significant at 0.01, GO tree depth < 10
            dfS = dfA.loc[(
                (dfA['p_fdr_bh'] <= 0.05) &
                #(dfA['depth'] < 10) &
                (dfA['NS'] == 'BP')),
            :]
            # Redo Index
            n = len(dfS)
            index = pd.MultiIndex.from_arrays([[celltype] * n, [layer] * n, [mid] * n, [mname] * n], names=('celltype', 'layer', 'module-id', 'module-name'))
            dfS.index = index

            ldfS.append(dfS)

    dfR = pd.concat(ldfS, axis=0)

    print("Exporting")
    ensurePathExists(wCSVFile)
    dfR.to_csv(wCSVFile)
