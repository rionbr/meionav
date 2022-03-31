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
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import precision_recall_fscore_support, matthews_corrcoef, roc_curve, precision_recall_curve, auc, confusion_matrix
from sklearn.dummy import DummyClassifier
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

    #
    # Feature Fix
    #

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
    cat_cols = ['biotype'] + ['mdlc-mutant-up/down']

    cat_features = df[cat_cols].fillna(0).T.to_dict().values()
    num_features = df[num_cols].values
    bin_features = df[bin_cols]

    DV = DictVectorizer(sparse=False)
    MMS = MinMaxScaler()
    #
    X_cat = DV.fit_transform(cat_features)
    X_num = MMS.fit_transform(num_features)
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

    classifiers = OrderedDict([
        ('Uniform Dummy', DummyClassifier(strategy='uniform', random_state=random_state)),
        ('Most Frequent Dummy', DummyClassifier(strategy='most_frequent', random_state=random_state)),
        ('Linear SVM', LinearSVC(max_iter=10000)),
        ('Random Forest', RandomForestClassifier()),
        ('Logistic Regression', LogisticRegression(solver='liblinear', max_iter=500)),
    ])

    # Classifier
    print('-- Fitting --')
    r = list()
    t = list()
    f = list()

    fold = 1
    for idx_train, idx_test in StratifiedKFold(n_splits=4, random_state=random_state, shuffle=True).split(X, y):
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

            precisions, recalls, thresholds = precision_recall_curve(y_test, probas_pred, pos_label=1)
            fprs, tprs, thresholds = roc_curve(y_test, probas_pred, pos_label=1)

            pr_auc = auc(recalls, precisions)
            roc_auc = auc(fprs, tprs)

            for y_test_, probas_pred_ in zip(y_test, probas_pred):
                t.append((clf_name, fold, y_test_, probas_pred_))

            # Eval
            precision, recall, f1, support = precision_recall_fscore_support(y_test, y_pred, average='binary', pos_label=1)
            mcc = matthews_corrcoef(y_test, y_pred)

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
    wMLTFile = 'results/clf-conserved-metric/{layer:s}/ml-{celltype:s}-{layer:s}-auc.csv.gz'.format(celltype=celltype, layer=layer)
    ensurePathExists(wMLTFile)
    dfT.to_csv(wMLTFile, encoding='utf-8', index=False)

    print('> Features')
    dfF = pd.DataFrame(f, columns=['clf', 'fold', 'feature', 'coef'])
    wMLFFile = 'results/clf-conserved-metric/{layer:s}/ml-{celltype:s}-{layer:s}-features.csv.gz'.format(celltype=celltype, layer=layer)
    ensurePathExists(wMLFFile)
    dfF.to_csv(wMLFFile, encoding='utf-8', index=True)

    print('> Classification')
    dfR = pd.DataFrame(r, columns=['clf', 'fold', 'TP', 'TN', 'FP', 'FN', 'precision', 'recall', 'f1', 'support', 'mcc', 'roc_auc', 'pr_auc'])
    wMLCFile = 'results/clf-conserved-metric/{layer:s}/ml-{celltype:s}-{layer:s}-clf.csv.gz'.format(celltype=celltype, layer=layer)
    ensurePathExists(wMLCFile)
    dfR.to_csv(wMLCFile, encoding='utf-8', index=False)
