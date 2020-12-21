# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
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
# from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.svm import LinearSVC #, SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
#from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
#from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction import DictVectorizer
#from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, cross_val_score
from sklearn.metrics import precision_recall_fscore_support, matthews_corrcoef, roc_curve, precision_recall_curve, auc, confusion_matrix
from sklearn.dummy import DummyClassifier
import gzip
from collections import OrderedDict
#
import argparse
from data import *


def get_module_dict_from_module_id(mid, modules):
    for module in modules:
        if module['id'] == mid:
            return module


def calc_y(r):
    # dfX.loc[dfX['core'] == True, ['core', 'mean-fert-rate', 'new-DM-phenotype', 'known-DM-phenotype', 'y']] 
    # Fertility Rate
    if not pd.isnull(r.get('mean-fert-rate', None)):
        if r.get('mean-fert-rate', 1.0) < 0.7:
            return True
        else:
            return False
    # Now Known
    # if not pd.isnull(r.get('new-DM-phenotype', None)):
    #   return True
    # Previously Known
    if not pd.isnull(r.get('known-DM-phenotype', None)):
        return True
    # Not tested yet
    if r.get('core', False):
        return None  # THIS IS THE SPECIAL CASE
    # ALL ELSE
    return False


if __name__ == '__main__':

    celltype = 'spermatocyte'
    layer = 'DM'
    mid = 1  # module id
    #
    modules = data_cells[celltype]['modules-svd']['modules'][layer]
    module = get_module_dict_from_module_id(mid, modules)
    name = module['name']
    #
    rMLFile = 'results/ml/{celltype:s}/{layer:s}/ml-{celltype:s}-{layer:s}-mod-{mid:d}.csv.gz'.format(celltype=celltype, layer=layer, mid=mid)
    dfX = pd.read_csv(rMLFile, index_col=0)
    #
    print(dfX.head())
    #
    # Feature Fix
    #
    #dfX.drop(labels=['degree', 'degree_weight', 'degree_backbone'], axis='columns', inplace=True)

    # Core
    c = 'core'
    if c in dfX.columns:
        dfX.loc[dfX[c].isnull(), c] = False

    # Validated-RNAi
    c = 'validated-rnai'
    if c in dfX.columns:
        dfX.loc[dfX[c].isin(['Yes']), c] = True
        dfX.loc[dfX[c].isin(['No', '-']), c] = False
        dfX.loc[dfX[c].isnull(), c] = False

    # Meiotic Entry
    c = 'meiotic-entry'
    if c in dfX.columns:
        dfX.loc[dfX[c].isnull(), c] = False

    # Meiotic Exit
    c = 'meiotic-exit'
    if c in dfX.columns:
        dfX.loc[dfX[c].isnull(), c] = False

    dfX['y'] = dfX.apply(calc_y, axis='columns')

    # Remove y == 'None' (those haven't been tested yet)
    dfX = dfX.loc[~dfX['y'].isnull() , :]

    #
    # Features
    #
    # Categorical
    if layer == 'DM':
        cat_HS_cols = []
        cat_MM_cols = []
        cat_DM_cols = [
            'spermatocyte',
            'biotype',
            'validated-rnai'
        ]
    elif layer in ['HS', 'MM']:
        cat_HS_cols = ['mammals']
        cat_MM_cols = ['mammals']
        cat_DM_cols = []
    else:
        raise Exception('Layer not defined')
    cat_cols = [] + cat_HS_cols + cat_MM_cols + cat_DM_cols
    # Numerical
    num_cols = [
        'TPM',
        'logFPKM',
        'degree',
        'degree_weight',
        'average_neighbor_degree',
        'degree_experiments_weight',
        'degree_backbone',
        'degree_centrality',
        'eigenvector_centrality',
        'betweenness_centrality',
        'pagerank',
        'clustering',
        'eccentricity',
        'core_number',
        'PCA-1c', 'PCA-2c', 'PCA-3c', 'PCA-4c', 'PCA-5c', 'PCA-6c', 'PCA-7c', 'PCA-8c', 'PCA-9c',
    ]
    bin_cols = [
        'meiotic-entry',
        'meiotic-exit'
    ]


    cat_features = dfX[cat_cols].T.to_dict().values()
    num_features = dfX[num_cols].values
    bin_features = dfX[bin_cols]


    DV = DictVectorizer(sparse=False)
    SS = StandardScaler()

    print('Vectorizing Categorical Features')
    X_cat = DV.fit_transform(cat_features)
    print('Scaling Numerical Features')
    X_num = SS.fit_transform(num_features)
    print('Binary Features')
    X_bin = bin_features

    features = DV.feature_names_ + num_cols + bin_cols

    print('Combine features')
    X = np.hstack((X_cat, X_num, X_bin))
    y = dfX['y'].astype(int).values

    n_instances, n_features = X.shape
    n_positive = len(y[y == True])
    n_negative = len(y[y == False])

    if n_positive <= 0:
        raise Exception("No positive features to predict. Sorry babe.")

    print('Instances: {:,d}'.format(n_instances))
    print('Features: {:,d}'.format(n_features))
    print()
    print('Positives: {:,d} ({:.2%})'.format(n_positive, (n_positive / n_instances)))
    print('Negatives: {:,d} ({:.2%})'.format(n_negative, (n_negative / n_instances)))
    print()

    random_state = 123

    classifiers = OrderedDict([
        ('Uniform Dummy', DummyClassifier(strategy='uniform', random_state=random_state)),
        ('Biased Dummy', DummyClassifier(strategy='stratified', random_state=random_state)),
        ('Linear SVM', LinearSVC(max_iter=10000)),
        ('Random Forest', RandomForestClassifier()),
        ('Logistic Regression', LogisticRegression(solver='liblinear', max_iter=500)),
    ])


    # Classifier
    print('-- Fitting --')
    r = list()
    t = list()
    f = list()

    random_state = 123
    fold = 1
    ##for idx_train, idx_test in StratifiedKFold(n_splits=2, random_state=random_state, shuffle=True).split(X, y):
    for idx_train, idx_test in LeaveOneOut().split(X, y):
        print('Split Fold: {fold:d}'.format(fold=fold))
        X_train, y_train, X_test, y_test = X[idx_train], y[idx_train], X[idx_test], y[idx_test]

        for clf_name, clf in classifiers.items():
            # print('Fitting: {clf_name:s}'.format(clf_name=clf_name))

            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)

            if hasattr(clf, 'predict_proba'):
                probas_pred = clf.predict_proba(X_test)[:, 1]
            else:
                probas_pred = clf.decision_function(X_test)
                # probas_pred = (probas_pred - probas_pred.min()) / (probas_pred.max() - probas_pred.min())

            """
            precisions, recalls, thresholds = precision_recall_curve(y_test, probas_pred, pos_label=1)
            fprs, tprs, thresholds = roc_curve(y_test, probas_pred, pos_label=1)

            pr_auc = auc(recalls, precisions)
            roc_auc = auc(fprs, tprs)
            """
            if 'Dummy' not in clf_name:
                # To CSV
                for y_test_, probas_pred_ in zip(y_test, probas_pred):
                    t.append((clf_name, fold, y_test_, probas_pred_))

            # Eval
            """
            precision, recall, f1, support = precision_recall_fscore_support(y_test, y_pred, average='binary', pos_label=1)
            mcc = matthews_corrcoef(y_test, y_pred)
            """
            TN, FP, FN, TP = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()

            # Features
            if hasattr(clf, 'coef_'):
                coefs = clf.coef_ if len(clf.coef_) > 1 else clf.coef_[0]
                for feature, coef in zip(features, coefs):
                    f.append((clf_name, fold, feature, coef))

            # Results
            r.append((clf_name, fold, TP, TN, FP, FN, precision, recall, f1, support, mcc, roc_auc, pr_auc))
        fold += 1

    ##
    # Export
    ##
    print('Export Results')

    print('> AUC Thresholds')
    dfT = pd.DataFrame(t, columns=['clf', 'fold', 'y_test', 'probas'])
    wMLTFile = 'results/ml/ml_thresholds_{layer:s}-SVD-{col:s}.csv.gz'.format(layer=layer, col=col)
    dfT.to_csv(wMLTFile, encoding='utf-8', index=False)

    print('> Features')
    dfF = pd.DataFrame(f, columns=['clf', 'fold', 'feature', 'coef'])
    dfF = dfF.groupby(['clf', 'feature']).agg({'coef': 'mean'})
    wMLFFile = 'results/ml/ml_features_{layer:s}-SVD-{col:s}.csv.gz'.format(layer=layer, col=col)
    dfF.to_csv(wMLFFile, encoding='utf-8', index=True)

    print('> Classification')
    dfR = pd.DataFrame(r, columns=['clf', 'fold', 'TP', 'TN', 'FP', 'FN', 'precision', 'recall', 'f1', 'support', 'mcc', 'roc_auc', 'pr_auc'])
    dfR.sort_values('clf', inplace=True)
    dfR.groupby('clf').agg('sum')
    wMLCFile = 'results/ml/ml_clfs_{layer:s}-SVD-{col:s}.csv.gz'.format(layer=layer, col=col)
    dfR.to_csv(wMLCFile, encoding='utf-8', index=False)
    print(dfR.groupby('clf').agg('mean'))
