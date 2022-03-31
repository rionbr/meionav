
# coding=utf-8
# Author: Rion B Correia
# Date: Sept 02, 2019
#
# Description: Reads a MultiLayer network (HS, MM & DM) and extracts subgraphs based on parameters for the networkbrowser.
#
#
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import networkx as nx
from matplotlib import colors
from utils import get_network_layer, get_network_largest_connected_component, ensurePathExists
import argparse
#
from data_spermatocyte_pca_modules_dm import spermatocyte_pca_modules_dm
from data_spermatocyte_pca_modules_mm import spermatocyte_pca_modules_mm
from data_spermatocyte_pca_modules_hs import spermatocyte_pca_modules_hs
#
from data_enterocyte_pca_modules_dm import enterocyte_pca_modules_dm
from data_enterocyte_pca_modules_mm import enterocyte_pca_modules_mm
from data_enterocyte_pca_modules_hs import enterocyte_pca_modules_hs


cmap_meanfertrate = colors.LinearSegmentedColormap.from_list(name='cmap-mean-fert-rate', colors=['#d62728', '#1f77b4'], N=256)


def fert_rate_color(x):
    if pd.isnull(x):
        return '#FFFFFF'  # white
    else:
        return colors.to_hex(cmap_meanfertrate(x))  # color


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=['spermatocyte', 'enterocyte'], help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    #
    parser.add_argument("--add_modules", default=True, type=bool, help="Add PCA module information to the network.")
    parser.add_argument("--add_conserved", default=True, type=bool, help="Add gene conservation information to the network.")
    parser.add_argument("--add_core", default=True, type=bool, help="Add core gene information to the network.")
    parser.add_argument("--add_backbone", default=True, type=bool, help="Add edge backbone to the network.")
    parser.add_argument("--add_ortho_backbone", default=True, type=bool, help="Add edge ortho-backbone to the network.")
    #
    parser.add_argument("--add_mdlc_dge_results", default=True, type=bool, help="Add gene mdlc DGE results to the DM network.")
    parser.add_argument("--add_splicing_defects", default=True, type=bool, help="Add gene mdlc splicing defects results to the DM network.")
    #
    parser.add_argument("--remove_isolates", default=True, type=bool, help="Remove isolate nodes from layers.")
    parser.add_argument("--only_largest_component", default=True, type=bool, help="Only output the largest connected component.")

    # parser.add_argument("--layer", default='DM', type=str, choices=['DM', 'MM', 'HS'], help="Network layer to compute SVD. Defaults to 'DM'.")

    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    #
    add_modules = args.add_modules
    add_conserved = args.add_conserved
    add_core = args.add_core
    add_backbone = args.add_backbone
    add_ortho_backbone = args.add_ortho_backbone
    #
    add_mdlc_dge_results = args.add_mdlc_dge_results
    add_splicing_defects = args.add_splicing_defects
    #
    remove_isolates = args.remove_isolates
    only_largest_component = args.only_largest_component
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
    print('Reading Network')
    path_net = '../../04-network/results/network/{celltype:}/'.format(celltype=celltype)
    if network == 'thr':
        rGfile_gpickle = path_net + 'net-{celltype:s}-{network:s}-{threshold:s}.gpickle'.format(celltype=celltype, network=network, threshold=threshold_str)
    G = nx.read_gpickle(rGfile_gpickle)

    if add_conserved:
        path_fpkm = '../../02-core_genes/results/FPKM/'
        df_HS = pd.read_csv(path_fpkm + 'HS/HS-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
        df_MM = pd.read_csv(path_fpkm + 'MM/MM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')
        df_DM = pd.read_csv(path_fpkm + 'DM/DM-FPKM-{celltype:s}.csv.gz'.format(celltype=celltype), index_col='id_string')

        dict_string_gene_HS = df_HS['id_gene'].to_dict()
        dict_string_gene_MM = df_MM['id_gene'].to_dict()
        dict_string_gene_DM = df_DM['id_gene'].to_dict()

        print('Loading {celltype:s} meta genes'.format(celltype=celltype))
        path = '../../02-core_genes/results/'
        dfM = pd.read_csv(path + 'meta-genes/meta-{celltype:s}-genes.csv.gz'.format(celltype=celltype), index_col='id_eggnog', usecols=['id_eggnog', 'id_string_HS', 'id_string_MM', 'id_string_DM'])

        dfM['id_string_HS'] = dfM['id_string_HS'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
        dfM['id_string_MM'] = dfM['id_string_MM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])
        dfM['id_string_DM'] = dfM['id_string_DM'].apply(lambda x: x.split(',') if not pd.isnull(x) else [])

        dfM['id_gene_HS'] = dfM['id_string_HS'].apply(lambda x: [dict_string_gene_HS[i] for i in x])
        dfM['id_gene_MM'] = dfM['id_string_MM'].apply(lambda x: [dict_string_gene_MM[i] for i in x])
        dfM['id_gene_DM'] = dfM['id_string_DM'].apply(lambda x: [dict_string_gene_DM[i] for i in x])

        dfM = dfM[['id_gene_HS', 'id_gene_MM', 'id_gene_DM']]
        # Only keep meta genes with homologs in all three species
        dfM = dfM.loc[dfM.applymap(len).applymap(bool).sum(axis='columns') == 3]

    for layer in ['HS', 'MM', 'DM']:
        if layer != 'DM':
            continue
        #
        print('Isolate {layer:s} Layer'.format(layer=layer))
        Gt = get_network_layer(G, layer=layer)

        # Add Module Information
        if add_modules:
            print('Load PCA Results ({layer:s})'.format(layer=layer))

            rPCAFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            rDiAnFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dian.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            rEntFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-entropy.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)

            df_pca = pd.read_csv(rPCAFile, index_col=0)
            df_dian = pd.read_csv(rDiAnFile, index_col=0)
            df_ent = pd.read_csv(rEntFile, index_col=0)

            modules = data[celltype]['modules'][layer]
            #
            rPCAFile = '../../04-network/results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            dfPCA = pd.read_csv(rPCAFile, index_col=0, encoding='utf-8')

            # Extract Component Modules
            for module in modules:

                mid = module['id']
                mname = module['name']
                print("Identifying module: M{mid:d} {mname:s} ".format(mid=mid, mname=mname))

                xc = module['dim-coords']['xdim']
                yc = module['dim-coords']['ydim']
                cx = "{xc:d}c".format(xc=xc)
                cy = "{yc:d}c".format(yc=yc)
                cxy = '{xc:d}c-{yc:d}c-dist'.format(xc=xc, yc=yc)  # label-1c-2c-dist
                cxl, cxh = module['dim-coords']['xvals']
                cyl, cyh = module['dim-coords']['yvals']
                cutrank = module['dim-coords']['radius-rank']
                cutradius = df_ent.loc[((df_ent['dim'] == xc) & (df_ent['cut-rank'] == cutrank)), 'radius-start'].squeeze()

                df_pca_tmp = df_pca.loc[
                    (
                        (df_pca[cx] >= cxl) & (df_pca[cx] <= cxh) & (df_pca[cy] >= cyl) & (df_pca[cy] <= cyh) & (df_dian[cxy] >= cutradius)
                    ), ['gene', cx, cy]].copy()

                component_ids = {g: True for g in df_pca_tmp.index.tolist()}
                net_attribute_name = 'module-pca-{layer:s}-{mid:d}'.format(layer=layer, mid=mid)
                nx.set_node_attributes(Gt, values=component_ids, name=net_attribute_name)

        # Add Conserved Information
        if add_conserved:
            dict_conserved = {gene: True for gene in dfM['id_gene_' + layer].explode().tolist()}
            #
            nx.set_node_attributes(Gt, values=dict_conserved, name='conserved')

        if add_core:
            rCOREFile = '../../02-core_genes/results/pipeline-core/{layer:s}_meiotic_genes.csv'.format(layer=layer)
            dfC = pd.read_csv(rCOREFile, index_col=0)
            dict_core = {gene: True for gene in dfC.index.tolist()}
            #
            nx.set_node_attributes(Gt, values=dict_core, name='core')

        # Remove Isolates
        if remove_isolates:
            isolates = list(nx.isolates(Gt))
            print('Removing {n:d} isolated nodes'.format(n=len(isolates)))
            Gt.remove_nodes_from(isolates)

        if add_backbone:
            print('Adding backbone data')
            path_backbone = "../../04-network/results/network-closure/{celltype:s}/".format(celltype=celltype)
            rBfile = path_backbone + "net-closure-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle".format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            B = nx.read_gpickle(rBfile)
            #
            is_metric = nx.get_edge_attributes(B, name='is_metric')
            is_ultrametric = nx.get_edge_attributes(B, name='is_ultrametric')
            #
            nx.set_edge_attributes(Gt, name='is_metric', values=is_metric)
            nx.set_edge_attributes(Gt, name='is_ultrametric', values=is_ultrametric)

        if add_ortho_backbone and (celltype == 'spermatocyte'):
            print('Adding ortho-backbone data')
            path_ortho_backbone = "../../04-network/results/network-closure-ortho/{celltype:s}/".format(celltype=celltype)
            rOBfile = path_ortho_backbone + "net-closure-ortho-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.gpickle".format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
            OB = nx.read_gpickle(rOBfile)
            #
            is_metric_ortho = nx.get_edge_attributes(OB, name='is_metric_ortho')
            #
            nx.set_edge_attributes(Gt, name='is_metric_ortho', values=is_metric_ortho)

        if add_mdlc_dge_results and (celltype == 'spermatocyte') and (layer == 'DM'):
            print('Adding mdlc DGE results')
            rMDLCFile = '../../01-diff-gene-exp/results/mdlc/{layer:s}-DGE-mdlc_vs_control.csv'.format(layer=layer)
            dfM = pd.read_csv(rMDLCFile, index_col=0, usecols=['id', 'gene', 'logFC', 'logCPM', 'F', 'PValue', 'FDR'])
            # Filter only DGE significant
            dfM = dfM.loc[(dfM['logFC'].abs() > 1) & (dfM['FDR'] <= 0.05) & (dfM['logCPM'] >= 1), :].copy()
            dfM['up/down'] = dfM['logFC'].map(lambda x: 'up' if x > 0 else 'down')
            #
            is_mdlc_dge = dfM['up/down'].to_dict()
            #
            nx.set_node_attributes(Gt, name='is_mdlc_dge', values=is_mdlc_dge)

        if add_splicing_defects and (celltype == 'spermatocyte') and (layer == 'DM'):
            print('Adding mdlc Splicing Defects results')
            rMDLCFile = '../../01-diff-gene-exp/results/mdlc/{layer:s}-IntronRetention-mdlc_vs_control.csv'.format(layer=layer)
            dfI = pd.read_csv(rMDLCFile, index_col=0, usecols=['id', 'gene'])
            #
            is_mdlc_intron = {n: True for n in dfI.index}
            #
            nx.set_node_attributes(Gt, name='is_mdlc_intron', values=is_mdlc_intron)


        # Largest Connected Component
        if only_largest_component:
            Gt = get_network_largest_connected_component(Gt)

        # graphml
        print('Export to graphml')
        if network == 'thr':
            wGtfile_graphml = '../gephi-plotting/results/graphml/{celltype:s}/net-{celltype:s}-{network:s}-{threshold:s}-{layer:s}.graphml'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
        ensurePathExists(wGtfile_graphml)
        nx.write_graphml(Gt, wGtfile_graphml)
