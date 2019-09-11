# coding=utf-8
# Author: Rion B Correia
# Date: Sept 11, 2019
#
# Description: Dispalys Statistics about the calculated pipelines.
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.precision', 4)
from tabulate import tabulate


def df2md(df, y_index=False, *args, **kwargs):
    blob = tabulate(df, headers='keys', tablefmt='pipe', *args, **kwargs)
    if not y_index:
        return '\n'.join(['| {}'.format(row.split('|', 2)[-1]) for row in blob.split('\n')])
    return blob

if __name__ == '__main__':

    #
    # All 3 Species Conserved FDRp05
    #
    pipeline = 'all3-conserved-FDRp05'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_DM = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

    n_m = df_M.shape[0]
    n_hs = df_HS.shape[0]
    n_mm = df_MM.shape[0]
    n_dm = df_DM.shape[0]

    df_stat = pd.DataFrame.from_records([
        ('Meta', n_m),
        ('HS', n_hs),
        ('MM', n_mm),
        ('DM', n_dm),
    ], columns=['Species', 'Genes'])
    print(df2md(df_stat, floatfmt='.4f'))
    print('\n')


    #
    # All 3 Species Conserved FDRp01
    #
    pipeline = 'all3-conserved-FDRp01'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_DM = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

    n_m = df_M.shape[0]
    n_hs = df_HS.shape[0]
    n_mm = df_MM.shape[0]
    n_dm = df_DM.shape[0]

    df_stat = pd.DataFrame.from_records([
        ('Meta', n_m),
        ('HS', n_hs),
        ('MM', n_mm),
        ('DM', n_dm),
    ], columns=['Species', 'Genes'])
    print(df2md(df_stat, floatfmt='.4f'))
    print('\n')

    #
    # All 3 Species Pooling FDRp05
    #
    pipeline = 'all3-pooling-DM-FDRp05'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_DM = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

    n_m = df_M.shape[0]
    n_hs = df_HS.shape[0]
    n_mm = df_MM.shape[0]
    n_dm = df_DM.shape[0]

    df_stat = pd.DataFrame.from_records([
        ('Meta', n_m),
        ('HS', n_hs),
        ('MM', n_mm),
        ('DM', n_dm),
    ], columns=['Species', 'Genes'])
    print(df2md(df_stat, floatfmt='.4f'))
    print('\n')

    #
    # All 3 Species Pooling FDRp01
    #
    pipeline = 'all3-pooling-DM-FDRp01'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_DM = pd.read_csv('results/{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

    n_m = df_M.shape[0]
    n_hs = df_HS.shape[0]
    n_mm = df_MM.shape[0]
    n_dm = df_DM.shape[0]

    df_stat = pd.DataFrame.from_records([
        ('Meta', n_m),
        ('HS', n_hs),
        ('MM', n_mm),
        ('DM', n_dm),
    ], columns=['Species', 'Genes'])
    print(df2md(df_stat, floatfmt='.4f'))
    print('\n')

    #
    # Mamals Conserved FDRp05
    #
    pipeline = 'mammals-conserved-FDRp05'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

    n_m = df_M.shape[0]
    n_hs = df_HS.shape[0]
    n_mm = df_MM.shape[0]

    df_stat = pd.DataFrame.from_records([
        ('Meta', n_m),
        ('HS', n_hs),
        ('MM', n_mm),
    ], columns=['Species', 'Genes'])
    print(df2md(df_stat, floatfmt='.4f'))
    print('\n')

    #
    # Mammals Pooling
    #
    pipeline = 'mammals-pooling'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

    n_m = df_M.shape[0]
    n_hs = df_HS.shape[0]
    n_MM = df_MM.shape[0]

    df_stat = pd.DataFrame.from_records([
        ('Meta', n_m),
        ('HS', n_hs),
        ('MM', n_mm),
    ], columns=['Species', 'Genes'])
    print(df2md(df_stat, floatfmt='.4f'))
    print('\n')