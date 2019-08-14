# coding=utf-8
# Author: Rion B Correia
# Date: Jul 24, 2019
#
# Description: Select DGE for DM
#
#
import math
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from pybiomart import Dataset
from utils import ensurePathExists


if __name__ == '__main__':

    #
    # [D]rosophila [M]elanogaster
    #

    # Query bioMart for Gene Name/Description
    DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
    dfDM = DM.query(attributes=['affy_drosophila_2', 'ensembl_gene_id'])
    dfDM = dfDM.dropna(subset=['AFFY Drosophila 2 probe']).set_index('AFFY Drosophila 2 probe')
    dfDM.index.name = 'Probe Set ID'

    print("> [D]rosophila [M]elanogaster")
    # Apical Testis (Up-Regulated)
    df = pd.read_csv('data/DM/DM_ApicalTestis_Exp_Pcall.csv', index_col=0)
    df = df.loc[(df['Combined detection signal'] == 'P'), ['Ave signal', 'Combined detection signal']]
    df.columns = ['avg-signal', 'A/P']
    df = pd.merge(df, dfDM, left_index=True, right_index=True, how='left').reset_index()
    df = df.dropna(subset=['Gene stable ID']).set_index('Gene stable ID')
    n = df.shape[0]
    print("Apical Testis: n={:,d} Up-Related Genes".format(n))
    wCSVFile = 'results/DM/DE/DM_DE_UpApicalTestis.csv'
    ensurePathExists(wCSVFile)
    df.to_csv(wCSVFile)

    # MidTestis (Down-Regulated)
    df = pd.read_csv('data/DM/DM_MidTestis_Exp_Pcall.csv', index_col=0)
    df = df.loc[(df['Combined detection signal'] == 'P'), ['Ave signal', 'Combined detection signal']]
    df.columns = ['avg-signal', 'A/P']
    df = pd.merge(df, dfDM, left_index=True, right_index=True, how='left').reset_index()
    df = df.dropna(subset=['Gene stable ID']).set_index('Gene stable ID')
    n = df.shape[0]
    print("Mid Testis: n={:,d} Down-Related Genes\n".format(n))
    wCSVFile = 'results/DM/DE/DM_DE_DownMidTestis.csv'
    ensurePathExists(wCSVFile)
    df.to_csv(wCSVFile)

    # Compute results for both FDR <= 0.05 and FDR <= 0.01
    for maxFDR in [0.01, 0.05]:
        print("\n- FDR <= {:f}\n".format(maxFDR))
        #
        # [H]omo [S]apiens
        #
        minLogFC = math.log2(2)
        #
        maxFDR_str = "FDR_{:s}".format(str(maxFDR).replace(".", "p"))

        print("> [H]omo [S]apiens")
        #
        # Cyte vs Gonia (Up-Regulated in Cytes)
        #
        df = pd.read_csv("results/HS/HS-DGE_Cyte_vs_Gonia.csv", index_col=0)
        df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
        dfU = df.loc[(df['logFC'] >= 0), :]

        n = dfU.shape[0]
        print("Cyte vs Gonia: n={:,d} Up-Related Genes in Cyte".format(n))
        wCSVFile = "results/HS/DE/HS_DE_UpCyte_vs_Gonia-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        ensurePathExists(wCSVFile)
        dfU.to_csv(wCSVFile)

        #
        # Cyte vs Tids (Down-Regulated in Tids)
        #
        df = pd.read_csv("results/HS/HS-DGE_Cyte_vs_Tid.csv", index_col=0)
        df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
        dfD = df.loc[(df['logFC'] <= 0), :]
        n = dfD.shape[0]
        print("Cyte vs Tid:   n={:,d} Down-Related Genes in Cyte\n".format(n))
        wCSVFile = "results/HS/DE/HS_DE_DownCyte_vs_Tid-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        dfD.to_csv(wCSVFile)

        print("> [M]us [M]usculus")
        #
        # Cyte vs Gonia (Up-Regulated in Cytes)
        #
        df = pd.read_csv("results/MM/MM-DGE_Cyte_vs_Gonia.csv", index_col=0)
        df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
        dfU = df.loc[(df['logFC'] >= 0), :]
        n = dfU.shape[0]
        print("Cyte vs Gonia: n={:,d} Up-Related Genes in Cyte".format(n))
        wCSVFile = "results/MM/DE/MM_DE_UpCyte_vs_Gonia-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        ensurePathExists(wCSVFile)
        dfU.to_csv(wCSVFile)

        #
        # Cyte vs Tids (Down-Regulated in Tids)
        #
        df = pd.read_csv("results/MM/MM-DGE_Cyte_vs_Tid.csv", index_col=0, nrows=None)
        df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
        dfD = df.loc[(df['logFC'] <= 0), :]
        n = dfD.shape[0]
        print("Cyte vs Tid:   n={:,d} Down-Related Genes in Cyte".format(n))
        wCSVFile = "results/MM/DE/MM_DE_DownCyte_vs_Tid-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        dfD.to_csv(wCSVFile)

    print("Done.")
