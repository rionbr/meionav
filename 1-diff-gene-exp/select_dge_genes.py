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

    minLogFC = math.log2(2)

    # Compute results for both FDR <= 0.05 and FDR <= 0.01
    for maxFDR in [0.01, 0.05]:
        print("\n- FDR <= {:f}\n".format(maxFDR))
        maxFDR_str = "FDR_{:s}".format(str(maxFDR).replace(".", "p"))
        #
        # [D]rosophila [M]elanogaster
        #
        print("> [D]rosophila [M]elanogaster")
        
        # Query bioMart for Gene Name/Description
        """
        DM = Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org')
        dfDM = DM.query(attributes=['affy_drosophila_2', 'ensembl_gene_id'])
        dfDM = dfDM.dropna(subset=['AFFY Drosophila 2 probe']).set_index('AFFY Drosophila 2 probe')
        dfDM.index.name = 'Probe Set ID'
        """
        #
        # Mid vs Apical Testis (Up-Regulated)
        #
        df = pd.read_csv('results/DM/DM-DGE_Middle_vs_Apical.csv', index_col=0)
        df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
        dfU = df.loc[(df['logFC'] >= 0), :]
        n = dfU.shape[0]
        print("Apical vs Mid Testis: n={:,d} Up-Related Genes".format(n))
        wCSVFile = "results/DM/DE/DM_DE_UpMiddle_vs_Apical-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        ensurePathExists(wCSVFile)
        dfU.to_csv(wCSVFile)

        #
        # MidTestis (Down-Regulated)
        #
        df = pd.read_csv('results/DM/DM-DGE_Middle_vs_Basal.csv', index_col=0)
        df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
        dfD = df.loc[(df['logFC'] <= 0), :]
        n = dfD.shape[0]
        print("Mid vs Distal Testis: n={:,d} Down-Related Genes\n".format(n))
        wCSVFile = "results/DM/DE/DM_DE_DownMiddle_vs_Basal-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        ensurePathExists(wCSVFile)
        dfD.to_csv(wCSVFile)


        #
        # [H]omo [S]apiens
        #
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
        print("Cyte vs Tid: n={:,d} Down-Related Genes in Cyte\n".format(n))
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
        print("Cyte vs Tid: n={:,d} Down-Related Genes in Cyte".format(n))
        wCSVFile = "results/MM/DE/MM_DE_DownCyte_vs_Tid-{maxFDR:s}.csv".format(maxFDR=maxFDR_str)
        dfD.to_csv(wCSVFile)

    print("Done.")
