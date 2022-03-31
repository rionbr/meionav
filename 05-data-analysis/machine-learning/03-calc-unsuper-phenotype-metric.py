# coding=utf-8
# Author: Rion B Correia
# Date: Jan 13, 2020
#
# Description: Reads a gene-feature table (X) and computed machine learning models
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
#
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_extraction import DictVectorizer
#
from utils import ensurePathExists
from collections import OrderedDict
import argparse


def get_module_dict_from_module_id(mid, modules):
    for module in modules:
        if module['id'] == mid:
            return module


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument('--layer', default='DM', type=str, choices=['HS', 'MM', 'DM'], help="Layer/Species.")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    layer = species = args.layer
    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')

    print('Load X-y features')
    rMLFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-X-y.csv.gz'.format(celltype=celltype, layer=layer)
    df = pd.read_csv(rMLFile, index_col=0)

    print('Load G features')
    rGFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-adj-matrix.csv.gz'.format(celltype=celltype, layer=layer)
    dfG = pd.read_csv(rGFile, index_col=0)

    print('Load B features')
    rBFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-adj-matrix-metric.csv.gz'.format(celltype=celltype, layer=layer)
    dfB = pd.read_csv(rBFile, index_col=0)
    #
    # Select Only core genes
    #
    dfG = dfG.loc[(df['core'].astype(bool).tolist()), :]
    dfB = dfB.loc[(df['core'].astype(bool).tolist()), :]

    # Only edge features with links
    dfG = dfG.loc[:, dfG.sum(axis='index').astype(bool)]
    dfB = dfB.loc[:, dfB.sum(axis='index').astype(bool)]

    df = df.loc[(df['core'].astype(bool)), :]
    #
    # Features
    #

    cols_logFC = ['logFC_MiddleApical', 'logFC_BasalMiddle'] if species == 'DM' else ['logFC_CyteGonia', 'logFC_TidCyte']
    cols_cellcomp = ['Middle_vs_Apical', 'Basal_vs_Middle'] if species == 'DM' else ['Cyte_vs_Gonia', 'Tid_vs_Cyte']

    if species == 'DM':
        cols_mods = ['module-1', 'module-2', 'module-3', 'module-4', 'module-5', 'module-6', 'module-11', 'module-12']
    elif species == 'MM':
        cols_mods = ['module-1', 'module-2', 'module-3', 'module-4', 'module-5', 'module-6', 'module-7', 'module-8']
    elif species == 'HS':
        cols_mods = ['module-1', 'module-2', 'module-3', 'module-4', 'module-5', 'module-7', 'module-8', 'module-10']
    #
    num_cols = [
        'TPM',
        'logFPKM',
        'degree',
        'degree-weight',
        'degree-centrality',
        'eigenvector-centrality',
        'pagerank',
        'avg-neighbor-degree'] + \
        cols_logFC + \
        [
        'logFC_mdlc-mutant',
        'degree-metric',
        'degree-weight-metric',
        'degree-centrality-metric',
        'clustering-metric',
        'pagerank-metric'
    ]
    bin_cols = cols_cellcomp + cols_mods + ['logFC_mdlc-mutant', 'mdlc-mutant-splidef']
    cat_cols = ['biotype']# + ['mdlc-mutant-up/down']

    cat_features = df[cat_cols].fillna(0).T.to_dict().values()
    num_features = df[num_cols].values
    amG_features = dfG.values
    amB_features = dfB.values
    bin_features = df[bin_cols]

    DV = DictVectorizer(sparse=False)
    MMS = MinMaxScaler()
    #
    X_cat = DV.fit_transform(cat_features)
    X_num = MMS.fit_transform(num_features)
    X_amg = MMS.fit_transform(amG_features)
    X_amb = MMS.fit_transform(amB_features)
    X_bin = bin_features
    #
    X = np.hstack((X_cat, X_num, X_bin))
    #y = df['phenotype'].fillna(0).astype(bool).values
    #
    n_instances, n_features = X.shape

    print('Instances: {:,d}'.format(n_instances))
    print('Features: {:,d}\n'.format(n_features))
    #
    random_state = 123

    # Classifier
    print('-- Fitting --')

    tsne = TSNE(n_components=2, early_exaggeration=12, learning_rate=200, n_iter=2000, init='pca', random_state=random_state)
    X_embed = tsne.fit_transform(X)

    ##
    # Export
    ##
    print('Export Results')

    print('> XY positions')
    dfXY = pd.DataFrame(X_embed, columns=['x', 'y'])
    wXYFile = 'results/unsuper-tsne/ml-unsuper-{celltype:s}-{layer:s}-tsne.csv.gz'.format(celltype=celltype, layer=layer)
    ensurePathExists(wXYFile)
    dfXY.to_csv(wXYFile, encoding='utf-8', index=False)
