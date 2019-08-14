# coding=utf-8
# Author: Rion B Correia
# Date: Jul 24, 2019
#
# Description: Select DGE for both HS and MM
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


if __name__ == '__main__':

    #
    # [H]omo [S]apiens
    #
    print("> [H]omo [S]apiens")
    #
    # Cyte vs Gonia (Up-Regulated in Cytes)
    #
    df = pd.read_csv("results/HS-DGE_Cyte_vs_Gonia.csv", index_col=0)
    df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
    dfU = df.loc[(df['logFC'] >= 0), :]

    n = dfU.shape[0]
    print("Cyte vs Gonia: n={:,d} Up-Related Genes in Cyte".format(n))
    dfU.to_csv("results/HS_DGE_UpCyte_vs_Gonia.csv")

    #
    # Cyte vs Tids (Down-Regulated in Tids)
    #
    df = pd.read_csv("results/HS-DGE_Cyte_vs_Tid.csv", index_col=0)
    df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
    dfD = df.loc[(df['logFC'] <= 0), :]
    n = dfD.shape[0]
    print("Cyte vs Tid:   n={:,d} Down-Related Genes in Cyte\n".format(n))
    dfD.to_csv("results/HS_DGE_DownCyte_vs_Tid.csv")

    print("> [M]us [M]usculus")
    #
    # Cyte vs Gonia (Up-Regulated in Cytes)
    #
    df = pd.read_csv("results/MM-DGE_Cyte_vs_Gonia.csv", index_col=0)
    df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
    dfU = df.loc[(df['logFC'] >= 0), :]
    n = dfU.shape[0]
    print("Cyte vs Gonia: n={:,d} Up-Related Genes in Cyte".format(n))
    dfU.to_csv("results/MM_DGE_UpCyte_vs_Gonia.csv")

    #
    # Cyte vs Tids (Down-Regulated in Tids)
    #
    df = pd.read_csv("results/MM-DGE_Cyte_vs_Tid.csv", index_col=0, nrows=None)
    df = df.loc[((df['FDR'] <= maxFDR) & (df['logFC'].abs() >= minLogFC)), :]
    dfD = df.loc[(df['logFC'] <= 0), :]
    dfG = pd.read_csv("goldstandard/MM_MeioticGenes_DownRegulated.csv", index_col=0)
    n = dfD.shape[0]
    print("Cyte vs Tid:   n={:,d} Down-Related Genes in Cyte".format(n))
    dfD.to_csv("results/MM_DGE_DownCyte_vs_Tid.csv")

    print("Done.")