# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Exports the list of genes belonging to the pca-modules identified
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
#
from data_spermatocyte_pca_modules_dm import spermatocyte_pca_modules_dm
from data_spermatocyte_pca_modules_mm import spermatocyte_pca_modules_mm
from data_spermatocyte_pca_modules_hs import spermatocyte_pca_modules_hs
#
from data_enterocyte_pca_modules_dm import enterocyte_pca_modules_dm
from data_enterocyte_pca_modules_mm import enterocyte_pca_modules_mm
from data_enterocyte_pca_modules_hs import enterocyte_pca_modules_hs


def export_genes(celltype='spermatocyte', network='thr', threshold=0.5, layer='DM', modules=[]):
    """ Export Genes """
    threshold_str = str(threshold).replace('.', 'p')
    #
    print('Exporting genes. PCA modules of {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
    rPCAFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rDiAnFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dian.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    rEntFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-entropy.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    #
    wCSVFile = 'results/pca-entropy/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-modules.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

    df_pca = pd.read_csv(rPCAFile, index_col=0)
    df_dian = pd.read_csv(rDiAnFile, index_col=0)
    df_ent = pd.read_csv(rEntFile, index_col=0)
    # df_cp = df_ent.loc[df_ent['cut-rank'].notnull(), :].sort_values(['dim', 'cut-rank'])

    ldfS = []
    for module in modules:

        mid = module['id']
        mname = module['name']
        print("Computing Module {mid:d}-{mname:s}".format(mid=mid, mname=mname))
        #
        xc = module['dim-coords']['xdim']
        yc = module['dim-coords']['ydim']
        ld1 = '{xc:d}c'.format(xc=xc)  # label 1 component
        ld2 = '{yc:d}c'.format(yc=yc)  # label 2 component
        l12d = '{xc:d}c-{yc:d}c-dist'.format(xc=xc, yc=yc)  # label-1c-2c-dist

        x0, x1 = module['dim-coords']['xvals']
        y0, y1 = module['dim-coords']['yvals']
        cut_rank = module['dim-coords']['radius-rank']

        # Radius of the circle
        cut_radius = df_ent.loc[((df_ent['dim'] == xc) & (df_ent['cut-rank'] == cut_rank)), 'radius-start'].squeeze()

        # Select genes in module
        df_pca_tmp = df_pca.loc[((df_pca[ld1] > x0) & (df_pca[ld1] < x1) & (df_pca[ld2] > y0) & (df_pca[ld2] < y1) & (df_dian[l12d] > cut_radius)), :].copy()
        #
        df_pca_tmp['module-id'] = mid
        df_pca_tmp['module-name'] = mname

        ldfS.append(df_pca_tmp)

    dfR = pd.concat(ldfS, axis=0)
    # Export
    print("Exporting")
    ensurePathExists(wCSVFile)
    dfR.to_csv(wCSVFile)


if __name__ == '__main__':

    celltype = 'enterocyte'  # spermatocyte or enterocyte
    network = 'thr'  # 'thr'
    threshold = 0.5
    layer = 'HS'

    data = {
        'spermatocyte': {
            'HS': spermatocyte_pca_modules_hs,
            'MM': spermatocyte_pca_modules_mm,
            'DM': spermatocyte_pca_modules_dm,
        },
        'enterocyte': {
            'HS': enterocyte_pca_modules_hs,
            'MM': enterocyte_pca_modules_mm,
            'DM': enterocyte_pca_modules_dm,
        }
    }
    modules = data[celltype][layer]

    export_genes(celltype, network, threshold, layer, modules)

    for celltype in ['spermatocyte', 'enterocyte']:
        for layer in ['HS', 'MM', 'DM']:
            modules = data[celltype][layer]
            export_genes(celltype, network, threshold, layer, modules)
