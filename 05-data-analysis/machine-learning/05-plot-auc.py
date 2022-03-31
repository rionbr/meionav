# coding=utf-8
# Author: Rion B Correia
# Date: Jan 14, 2021
#
# Description: Plot the AUC curves for the ML results
#
# Instructions:
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
#
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.handlelength'] = 1.5
#
mpl.rc('font', size=10)  # controls default text sizes
mpl.rc('axes', titlesize=16)  # fontsize of the axes title
mpl.rc('axes', labelsize=12)  # fontsize of the x and y labels
mpl.rc('xtick', labelsize=10)  # fontsize of the tick labels
mpl.rc('ytick', labelsize=10)  # fontsize of the tick labels
mpl.rc('legend', fontsize=10)  # legend fontsize
mpl.rc('figure', titlesize=16)  # fontsize of the figure title
#
mpl.rc('figure.subplot', left=0.15, right=0.9, bottom=0.12, top=0.89)
#
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, precision_recall_curve, auc
from sklearn.metrics import plot_precision_recall_curve, average_precision_score

from itertools import cycle
from utils import ensurePathExists


def plot_pr_auc(celltype, layer, clftype):
    print("-- Plot P/R AUC --")

    rCSVfile = 'results/clf-{clftype:s}/{layer:s}/ml-{celltype:s}-{layer:s}-auc.csv.gz'.format(clftype=clftype, celltype=celltype, layer=layer)
    df = pd.read_csv(rCSVfile)
    i_clf = 1
    i_clf = 1

    fig, ax = plt.subplots(figsize=(3.9, 3.7))

    #ax.set(adjustable='box', aspect='equal')  # equal weight and height

    clf_args = {
        'Random Forest': {'color': '#9467bd', 'linestyle': '-', 'zorder': 2},
        #'Logistic Regression': {'color': '#1f77b4', 'linestyle': '-', 'zorder': 4},
        'Linear SVM': {'color': '#e377c2', 'linestyle': '-', 'zorder': 3},
        #
        'Uniform Dummy': {'name': 'Coin toss', 'color': '#7f7f7f', 'linestyle': '--', 'zorder': 3, 'alpha': 0.7},
        #'Most Frequent Dummy': {'name': 'Majority vote', 'color': '#17becf', 'linestyle': (1, (3, 2)), 'zorder': 3, 'alpha': 0.7},
        #
        #'Baseline Dummy': {'name': 'Baseline', 'color': '#7f7f7f', 'linestyle': '--', 'zorder': 3}
    }
    clfs = clf_args.keys()
    for clf in clfs:
        print('-- Clf: {clf:s} --'.format(clf=clf))
        #
        dfc = df.loc[(df['clf'] == clf), :]
        args = clf_args.get(clf)
        #
        if args.get('name', None) is not None:
            clf_str = args.pop('name', None)
        else:
            clf_str = clf

        i_fold = 1

        mean_recall = np.linspace(0, 1, 100)
        precisions = []
        auc_prs = []

        for fold, dfcf in dfc.groupby('fold'):

            y_test = dfcf['y_test']
            probas_pred = dfcf['probas']

            pos_neg = y_test.value_counts()
            # n_positive = pos_neg.loc[1]
            # n_negative = pos_neg.loc[0]

            precision, recall, thresholds = precision_recall_curve(y_test, probas_pred, pos_label=1)
            avg_pr_score = average_precision_score(y_test, probas_pred, pos_label=1)
            # pr_auc = auc(recall, precision)

            precision = precision[::-1]
            recall = recall[::-1]
            thresholds = thresholds[::-1]

            if 'Dummy' in clf:
                precision[0] = precision[1]

            interp_precisions = np.interp(mean_recall, recall, precision)


            if 'Dummy' not in clf:
                interp_precisions[0] = 1
            precisions.append(interp_precisions)
            auc_prs.append(avg_pr_score)
            #
            i_fold += 1

        # Mean P/C
        mean_precision = np.mean(precisions, axis=0)
        mean_avg_pr_score = np.mean(auc_prs)
        # std_pr_auc = np.std(auc_prs)
        ax.plot(mean_recall, mean_precision, label='{clf:s} ({area:.2f})'.format(clf=clf_str, area=mean_avg_pr_score), lw=2, **args)

        # Fill
        if 'Dummy' not in clf:
            std_precision = np.std(precisions, axis=0)
            precisions_upper = np.minimum(mean_precision + std_precision, 1)
            precisions_lower = np.maximum(mean_precision - std_precision, 0)
            ax.fill_between(mean_recall, precisions_lower, precisions_upper, color=args['color'], alpha=.35, zorder=2)

        i_clf += 1

    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])

    ax.set_title('Precision & Recall Curve')
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    #
    ax.legend(loc="lower left")

    #plt.subplots_adjust(left=0.15, bottom=0.09, right=0.999, top=0.94, wspace=0.25, hspace=0.35)

    wIMGfile = 'images/ml-auc-{clftype:s}/ml-{celltype:s}-{clftype:s}-{layer:s}-pr-auc.pdf'.format(clftype=clftype, celltype=celltype, layer=layer)
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=300, pad_inches=0.0)


def plot_roc_auc(celltype, layer, clftype):
    print('-- Plot ROC AUC --')

    rCSVfile = 'results/clf-{clftype:s}/{layer:s}/ml-{celltype:s}-{layer:s}-auc.csv.gz'.format(clftype=clftype, celltype=celltype, layer=layer)
    df = pd.read_csv(rCSVfile)
    i_clf = 1

    fig, ax = plt.subplots(figsize=(3.9, 3.7))

    #ax.set(adjustable='box', aspect='equal')  # equal weight and height

    clf_args = {
        'Random Forest': {'color': '#9467bd', 'linestyle': '-', 'zorder': 2},
        #'Logistic Regression': {'color': '#1f77b4', 'linestyle': '-', 'zorder': 4},
        'Linear SVM': {'color': '#e377c2', 'linestyle': '-', 'zorder': 3},
        #
        'Uniform Dummy': {'name': 'Coin toss', 'color': '#7f7f7f', 'linestyle': '--', 'zorder': 3, 'alpha': 0.7},
        #'Most Frequent Dummy': {'name': 'Majority vote', 'color': '#17becf', 'linestyle': (1, (3, 2)), 'zorder': 3, 'alpha': 0.7},
        #
        #'Baseline Dummy': {'name': 'Baseline', 'color': '#7f7f7f', 'linestyle': '--', 'zorder': 3}
    }
    clfs = clf_args.keys()

    for clf in clfs:
        #
        print('-- Clf: {clf:s} --'.format(clf=clf))
        #
        dfc = df.loc[(df['clf'] == clf), :]
        args = clf_args.get(clf)
        #
        if args.get('name', None) is not None:
            clf_str = args.pop('name', None)
        else:
            clf_str = clf

        i_fold = 1

        mean_fpr = np.linspace(0, 1, 100)
        tprs = []
        auc_rocs = []
        #precisions = []

        for fold, dfcf in dfc.groupby('fold'):

            y_test = dfcf['y_test']
            probas_pred = dfcf['probas']
            #
            fpr, tpr, thresholds = roc_curve(y_test, probas_pred, pos_label=1)
            #
            roc_auc = auc(fpr, tpr)
            #
            #interp_precisions = np.interp(mean_recall, recall, precision)
            
            #if 'Dummy' not in clf:
            #    interp_precisions[0] = 1
            
            #precisions.append(interp_precisions)
            #
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            auc_rocs.append(roc_auc)

            i_fold += 1

        # Mean ROC
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_roc_auc = auc(mean_fpr, mean_tpr)
        # std_roc_auc = np.std(auc_rocs)
        ax.plot(mean_fpr, mean_tpr, label='{clf:s} ({area:.2f})'.format(clf=clf_str, area=mean_roc_auc), lw=2, **args)

        # Fill
        if 'Dummy' not in clf:
            std_tpr = np.std(tprs, axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color=args['color'], alpha=.35, zorder=2)

        i_clf += 1

    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    #
    ax.set_title('ROC Curve')
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')

    ax.legend(loc="lower right")

    #plt.subplots_adjust(left=0.15, bottom=0.09, right=0.999, top=0.94, wspace=0.25, hspace=0.35)

    wIMGfile = 'images/ml-auc-{clftype:s}/ml-{celltype:s}-{clftype:s}-{layer:s}-roc-auc.pdf'.format(clftype=clftype, celltype=celltype, layer=layer)
    ensurePathExists(wIMGfile)
    plt.savefig(wIMGfile, dpi=300, pad_inches=0.0)


if __name__ == '__main__':

    celltype = 'spermatocyte'
    layer = 'DM'
    clftype = 'phenotype'  # ['conserved', 'phenotype', 'conserved-without-exp']

    for layer in ['HS', 'MM', 'DM']:
        plot_pr_auc(celltype, layer, clftype)
        plot_roc_auc(celltype, layer, clftype)
