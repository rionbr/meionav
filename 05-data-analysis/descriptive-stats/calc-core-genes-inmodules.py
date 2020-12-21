# coding=utf-8
# Author: Rion B Correia
# Date: Dec 04, 2020
#
# Description: Loads modules and core genes and calculates stats
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists




if __name__ == '__main__':

    celltype = 'spermatocyte'
    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    layer = species = 'DM'

    #
    r = []
    for species in ['DM', 'MM', 'HS']:
        layer = species

        rCSVmFile = '../entropy-based-modules/results/pca-entropy/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-modules.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        rCSVcFile = '../../02-core_genes/results/pipeline-core/{species:s}_meiotic_genes.csv'.format(species=species)
        #
        dfm = pd.read_csv(rCSVmFile, index_col=0)
        dfc = pd.read_csv(rCSVcFile, index_col=0)
        cores = dfc.index.tolist()

        for (mid, mname), dfmt in dfm.groupby(['module-id', 'module-name']):
            ngenes = len(dfmt)
            ncoregenes = len(dfmt.loc[dfmt.index.isin(dfc.index)])

            r.append((celltype, layer, mid, mname, ngenes, ncoregenes))


    dfR = pd.DataFrame(r, columns=['celltype', 'species', 'module-id', 'module-name', 'n-genes', 'n-core-genes'])

    # Export
    dfR.to_csv('results/stats-core-genes-in-modules.csv')