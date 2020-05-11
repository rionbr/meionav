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
    # Mamals
    #
    pipeline = 'mammals'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/pipeline-{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/pipeline-{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/pipeline-{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

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
    # Core
    #
    pipeline = 'core'
    print('# Pipeline: {pipeline:s}\n'.format(pipeline=pipeline))

    df_M = pd.read_csv('results/pipeline-{pipeline:s}/meta_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_HS = pd.read_csv('results/pipeline-{pipeline:s}/HS_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_MM = pd.read_csv('results/pipeline-{pipeline:s}/MM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)
    df_DM = pd.read_csv('results/pipeline-{pipeline:s}/DM_meiotic_genes.csv'.format(pipeline=pipeline), index_col=0)

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
