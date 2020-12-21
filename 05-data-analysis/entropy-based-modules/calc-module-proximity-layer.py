# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads the SVD results and calculats a similarity between SVD Modules across layers.
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


if __name__ == '__main__':

    celltype_i = 'enterocyte'
    celltype_j = 'spermatocyte'
    #
    network = 'thr'  # 'thr'
    threshold = 0.5
    threshold_str = str(threshold).replace('.', 'p')
    #
    layers = ['HS', 'MM', 'DM']
    #
    placeholder = {'HS': None, 'MM': None, 'DM': None}
    data = {
        'spermatocyte': {
            'PCA': dict(placeholder),
            'distance-angle': dict(placeholder),
            'entropy': dict(placeholder),
            'modules': {
                'HS': spermatocyte_pca_modules_hs,
                'MM': spermatocyte_pca_modules_mm,
                'DM': spermatocyte_pca_modules_dm,
            },
        },
        'enterocyte': {
            'PCA': dict(placeholder),
            'distance-angle': dict(placeholder),
            'entropy': dict(placeholder),
            'modules': {
                'HS': enterocyte_pca_modules_hs,
                'MM': enterocyte_pca_modules_mm,
                'DM': enterocyte_pca_modules_dm,
            }
        }
    }

    #
    for celltype in [celltype_i, celltype_j]:
        for layer in layers:

            print('Loading PCA/Entropy for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
            rPCAFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            rDiAnFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dian.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            rEntFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-entropy.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

            df_pca = pd.read_csv(rPCAFile, index_col=0)
            df_dian = pd.read_csv(rDiAnFile, index_col=0)
            df_ent = pd.read_csv(rEntFile, index_col=0)
            #
            data[celltype]['PCA'][layer] = df_pca
            data[celltype]['distance-angle'][layer] = df_dian
            data[celltype]['entropy'][layer] = df_ent

    r = []
    for layer in layers:
        print("Comparing {celltype_i:s} with {celltype_j:s} of layer: {layer:s}".format(celltype_i=celltype_i, celltype_j=celltype_j, layer=layer))

        df_pca_i = data[celltype_i]['PCA'][layer]
        df_pca_j = data[celltype_j]['PCA'][layer]
        df_dian_i = data[celltype_i]['distance-angle'][layer]
        df_dian_j = data[celltype_j]['distance-angle'][layer]
        df_ent_i = data[celltype_i]['entropy'][layer]
        df_ent_j = data[celltype_j]['entropy'][layer]
        modules_i = data[celltype_i]['modules'][layer]
        modules_j = data[celltype_j]['modules'][layer]

        for module_i, module_j in [(a, b) for a in modules_i for b in modules_j]:

            id_i = str(module_i['id'])
            id_j = str(module_j['id'])
            name_i = module_i['name']
            name_j = module_j['name']
            print("Comparing: {layer:s}-{id_i:s}-{name_i:s} with {layer:s}-{id_j:s}-{name_j:s}".format(
                id_i=id_i, id_j=id_j, layer=layer, name_i=name_i, name_j=name_j)
            )

            ixc = module_i['dim-coords']['xdim']
            iyc = module_i['dim-coords']['ydim']
            icx = "{xc:d}c".format(xc=ixc)
            icy = "{yc:d}c".format(yc=iyc)
            icxy = '{xc:d}c-{yc:d}c-dist'.format(xc=ixc, yc=iyc)  # label-1c-2c-dist
            icxl, icxh = module_i['dim-coords']['xvals']
            icyl, icyh = module_i['dim-coords']['yvals']
            icutrank = module_i['dim-coords']['radius-rank']
            icutradius = df_ent_i.loc[((df_ent_i['dim'] == ixc) & (df_ent_i['cut-rank'] == icutrank)), 'radius-start'].squeeze()

            df_pca_i_tmp = df_pca_i.loc[
                (
                    (df_pca_i[icx] >= icxl) & (df_pca_i[icx] <= icxh) & (df_pca_i[icy] >= icyl) & (df_pca_i[icy] <= icyh) & (df_dian_i[icxy] >= icutradius)
                ), ['gene', icx, icy]].copy()

            jxc = module_j['dim-coords']['xdim']
            jyc = module_j['dim-coords']['ydim']
            jcx = "{xc:d}c".format(xc=jxc)
            jcy = "{yc:d}c".format(yc=jyc)
            jcxy = '{xc:d}c-{y:d}c-dist'.format(xc=jxc, y=jyc)  # label-1c-2c-dist
            jcxl, jcxh = module_j['dim-coords']['xvals']
            jcyl, jcyh = module_j['dim-coords']['yvals']
            jcutrank = module_i['dim-coords']['radius-rank']
            jcutradius = df_ent_j.loc[((df_ent_j['dim'] == jxc) & (df_ent_j['cut-rank'] == jcutrank)), 'radius-start'].squeeze()

            df_pca_j_tmp = df_pca_j.loc[
                (
                    (df_pca_j[jcx] >= jcxl) & (df_pca_j[jcx] <= jcxh) & (df_pca_j[jcy] >= jcyl) & (df_pca_j[jcy] <= jcyh) & (df_dian_j[jcxy] >= jcutradius)
                ), ['gene', jcx, jcy]].copy()

            genes_i = df_pca_i_tmp.index.to_list()
            genes_j = df_pca_j_tmp.index.to_list()

            # Jaccard Proximity
            a = set(genes_i)
            b = set(genes_j)
            a_union_b = a.union(b)
            a_inter_b = a.intersection(b)
            dist = len(a_inter_b) / len(a_union_b)
            r.append((layer, celltype_i, celltype_j, id_i, id_j, name_i, name_j, dist))

    dfR = pd.DataFrame(r, columns=['layer', 'celltype-i', 'celltype-j', 'id-i', 'id-j', 'name-i', 'name-j', 'proximity'])

    ##
    # Export
    ##
    print('Exporting')
    wCSVfile = 'results/module-proximity/module-proximity-layers-{network:s}-{threshold:s}.csv.gz'.format(network=network, threshold=threshold_str)
    ensurePathExists(wCSVfile)
    dfR.to_csv(wCSVfile)
