# coding=utf-8
# Author: Rion B Correia
# Date: Jun 30, 2019
#
# Description: Indexes certain genes and exports their list.
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import argparse
from utils import ensurePathExists
from scipy.stats import ks_2samp
from itertools import combinations
import swifter


#  Separating by At Least One Match
def select_by_at_least_one_match(ilist, keeplist):
    # Only keep genes that are found in any of our gene list
    return [i for i in ilist if i in keeplist]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Defaults to spermatocyte")
    parser.add_argument("--biotype", default='protein_coding', type=str, choices=['protein_coding'], help="Filter nodes by biotype (e.g., protein-coding)")
    parser.add_argument("--attribute", default='TPM', type=str, help="Which attribute to plot. Defaults to 'TPM'.")
    # parser.add_argument("--log", default=True, type=bool, help="Transforms attribute into log2(attribute).")
    parser.add_argument("--minTPM", default=1, type=int, help="minLogTPM = math.log2(x). Defaults to 1.")
    args = parser.parse_args()

    celltype = args.celltype  # spermatocyte or enterocyte
    biotype = args.biotype
    attribute = args.attribute
    # log = args.log
    minTPM = args.minTPM

    print('Exporint {celltype:s}-{biotype:s}-{attribute:s}'.format(celltype=celltype, biotype=biotype, attribute=attribute))

    print('Loading {celltype:s} Files'.format(celltype=celltype))
    path = '../../02-core_genes/results/'
    df_HS = pd.read_csv(path + 'FPKM/HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')
    df_MM = pd.read_csv(path + 'FPKM/MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')
    df_DM = pd.read_csv(path + 'FPKM/DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_gene')

    # Remove Duplicates
    df_HS = df_HS.loc[~df_HS.index.duplicated(keep='first'), :]
    df_MM = df_MM.loc[~df_MM.index.duplicated(keep='first'), :]
    df_DM = df_DM.loc[~df_DM.index.duplicated(keep='first'), :]

    # minTPM
    if minTPM:
        df_HS = df_HS.loc[(df_HS['TPM'] >= minTPM), :]
        df_MM = df_MM.loc[(df_MM['TPM'] >= minTPM), :]
        df_DM = df_DM.loc[(df_DM['TPM'] >= minTPM), :]

    # Meta Genes
    print('Loading {celltype:s} meta genes'.format(celltype=celltype))
    dfM = pd.read_csv(path + 'meta-genes/meta-{celltype:s}-genes.csv.gz'.format(celltype=celltype), index_col='id_eggnog', usecols=['id_eggnog', 'id_string_HS', 'id_string_MM', 'id_string_DM'])

    dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
    dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

    # Remove meta-genes not in selected species genes (df_HS) ..
    print("Subselecting meta-genes")
    dfM['id_string_HS'] = dfM['id_string_HS'].swifter.apply(select_by_at_least_one_match, args=(df_HS['id_string'].to_list(), ))
    dfM['id_string_MM'] = dfM['id_string_MM'].swifter.apply(select_by_at_least_one_match, args=(df_MM['id_string'].to_list(), ))
    dfM['id_string_DM'] = dfM['id_string_DM'].swifter.apply(select_by_at_least_one_match, args=(df_DM['id_string'].to_list(), ))

    # Only keep meta genes with homologs in all three species
    dfM = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 3]

    # Set of homologs
    set_id_string_HS = set(dfM['id_string_HS'].explode().to_list())
    set_id_string_MM = set(dfM['id_string_MM'].explode().to_list())
    set_id_string_DM = set(dfM['id_string_DM'].explode().to_list())

    # Index by homologs
    df_HS = df_HS.loc[df_HS['id_string'].isin(set_id_string_HS), :]
    df_MM = df_MM.loc[df_MM['id_string'].isin(set_id_string_MM), :]
    df_DM = df_DM.loc[df_DM['id_string'].isin(set_id_string_DM), :]

    print("Exporting")
    for specie, dft in zip(['HS', 'MM', 'DM'], [df_HS, df_MM, df_DM]):
        wCSVfile = 'results/gene-fpkm-export-conserved/{celltype:s}/csv-{celltype:s}-conserved-{specie:s}.csv.gz'.format(celltype=celltype, specie=specie)
        ensurePathExists(wCSVfile)
        dft.to_csv(wCSVfile)
