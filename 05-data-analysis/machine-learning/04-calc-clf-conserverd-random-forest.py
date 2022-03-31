# coding=utf-8
# Author: Rion B Correia
# Date: Jan 13, 2020
#
# Description: Reads a gene-feature table (X) and computed machine learning models
#
#
import os
import subprocess
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
#
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_extraction import DictVectorizer
#
import argparse
from dtreeviz.trees import dtreeviz


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument('--layer', default='HS', type=str, choices=['HS', 'MM', 'DM'], help="Layer/Species.")
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    layer = species = args.layer
    clftype = 'conserved'

    data = {
        'HS': {
            'color': '#2ca02c'
        },
        'MM': {
            'color': '#7f7f7f'
        },
        'DM': {
            'color': '#ff7f0e'
        }
    }

    print('Load X-y features')
    rMLFile = 'results/matrix-x-vector-y/ml-{celltype:s}-{layer:s}-X-y.csv.gz'.format(clftype=clftype, celltype=celltype, layer=layer)
    df = pd.read_csv(rMLFile, index_col=0)

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
        'avg-neighbor-degree'] + cols_logFC
    bin_cols = cols_cellcomp + cols_mods
    cat_cols = ['biotype']

    cat_features = df[cat_cols].fillna(0).T.to_dict().values()
    num_features = df[num_cols].values
    bin_features = df[bin_cols]

    DV = DictVectorizer(sparse=False)
    #
    X_cat = DV.fit_transform(cat_features)
    X_num = num_features
    X_bin = bin_features
    #
    features = DV.feature_names_ + num_cols + bin_cols
    #
    X = np.hstack((X_cat, X_num, X_bin))
    y = df['conserved'].fillna(0).astype(bool).values
    #
    n_instances, n_features = X.shape
    n_positive = len(y[y == True])
    n_negative = len(y[y == False])

    if n_positive <= 0:
        raise Exception("No positive features to predict. Sorry babe.")

    print('Instances: {:,d}'.format(n_instances))
    print('Features: {:,d}\n'.format(n_features))
    #
    print('Positives: {:,d} ({:.2%})'.format(n_positive, (n_positive / n_instances)))
    print('Negatives: {:,d} ({:.2%})\n'.format(n_negative, (n_negative / n_instances)))

    random_state = 123

    # Random Forest
    clf = DecisionTreeClassifier(criterion='gini', max_depth=3, min_samples_split=50, max_features='auto')

    # Fit
    print('-- Fitting --')
    clf.fit(X, y)

    #
    color = data[species]['color']
    #
    print('-- Visualizing Tree --')
    viz = dtreeviz(
        clf, X, y,
        target_name="Conserved",
        feature_names=features,
        class_names=['False', 'True'],
        colors={
            'classes': [None, None, [color, '#d62728']],
        },
        orientation='TD')

    # Export
    print('-- Exporting --')
    print('> Export SVG')
    wSVGfile = 'images/ml-tree-{clftype:s}/img-decision-tree-{species:s}.svg'.format(clftype=clftype, species=species)
    viz.save(wSVGfile)

    # convert to pdf
    print('> Convert to PDF')
    wPDFfile = 'images/ml-tree-{clftype:s}/img-decision-tree-{species:s}.pdf'.format(clftype=clftype, species=species)
    path = os.path.abspath(os.getcwd())
    input = path + '/' + wSVGfile
    output = path + '/' + wPDFfile
    cmd = "/Applications/Inkscape.app/Contents/Resources/script {input:s} --without-gui --export-pdf={output:s}".format(input=input, output=output)
    subprocess.call(cmd, shell=True)

    print('Done.')