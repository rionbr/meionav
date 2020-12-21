# coding=utf-8
# Author: Rion B Correia
# Date: June 17, 2020
#
# Description: Calculates entropy-based on network PCA
#
#
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from cycler import cycler
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
from utils import ensurePathExists
from scipy import stats
import argparse


def in_ranges(x, bins):
    return [((x >= lower) & (x <= upper)) for lower, upper in bins]


def compute_entropy(df_dec,
                    radius_window=1.0,
                    radius_overlap=0.1,
                    angle_window=30,
                    angle_overlap=15,
                    min_points=10,
                    n_cut_points=3,
                    components=9):
    """ """
    df_dec = df_dec.copy()
    #
    angle_items = int(angle_window / angle_overlap)
    a = np.arange(-180, (181), angle_overlap)
    angle_bins = [(i, j) for i, j in zip(a[0:-angle_items], a[angle_items:])]
    n_bins = len(angle_bins)
    max_entropy = stats.entropy((np.ones(shape=n_bins) / n_bins), base=2)
    list_df_ent = []
    #
    for dim in range(1, (components + 1)):
        print('Computing projection: {dim1:d} vs {dim2:d}'.format(dim1=dim, dim2=(dim + 1)))
        #
        cx = str(dim) + 'c'
        cy = str(dim + 1) + 'c'
        dist_label = '{cx:s}-{cy:s}-dist'.format(cx=cx, cy=cy)
        angle_label = '{cx:s}-{cy:s}-angle'.format(cx=cx, cy=cy)
        df_dec[dist_label] = np.hypot(df_dec[cx], df_dec[cy])
        df_dec[angle_label] = np.degrees(np.arctan2(df_dec[cy], df_dec[cx]))
        #
        df_dec.sort_values(dist_label, ascending=True, inplace=True)
        radius_max = df_dec[dist_label].max()
        #
        radius_items = int(radius_window / radius_overlap)
        #
        b = np.arange(0, (radius_max + radius_overlap), radius_overlap)
        radius_intervals = [(s, e) for s, e in zip(b[0:-radius_items], b[radius_items:])]

        # Loop radius intervals
        r = []
        for radius_start, radius_end in radius_intervals:

            df_dian_tmp = df_dec.loc[(df_dec[dist_label] >= radius_start) & (df_dec[dist_label] <= radius_end), :]

            dfc = df_dian_tmp[angle_label].apply(lambda x: pd.Series(in_ranges(x, angle_bins), angle_bins))

            if len(dfc) > min_points:
                dfp = (dfc.sum(axis=0) / dfc.sum().sum()).rename('prob').to_frame()
                dfp['log2'] = dfp['prob'].apply(np.log2)
                #
                entropy = stats.entropy(dfp['prob'], base=2)

            else:
                entropy = np.nan

            entropy_norm = entropy / max_entropy
            r.append((dim, radius_start, radius_end, entropy, entropy_norm))

        #
        df_ent_tmp = pd.DataFrame(r, columns=['dim', 'radius-start', 'radius-end', 'entropy', 'entropy-norm'])

        # Interpolation
        df_ent_tmp['entropy-smooth'] = df_ent_tmp['entropy-norm'].interpolate(method='linear', limit_direction='both')
        # Rank
        df_ent_tmp['radius-rank'] = df_ent_tmp['radius-start'].rank(method='min')
        df_ent_tmp['entropy-rank'] = df_ent_tmp['entropy-norm'].rank(method='min')
        # Rank Sum
        df_ent_tmp['rank-sum'] = ((df_ent_tmp['radius-rank']) + (df_ent_tmp['entropy-rank']))

        # Define Cut Pointns
        cut_points = []
        # Index % Sort
        df_cp = df_ent_tmp.sort_values('rank-sum').loc[(df_ent_tmp['radius-start'] > 1.0), :]
        possible_rank = 1
        for possible_id, row in df_cp.iterrows():
            possible_value = row['radius-start']
            if not any([True if abs(possible_value - existing_value) <= 1.0 else False for existing_id, existing_value, existing_rank in cut_points]):
                cut_points.append((possible_id, possible_value, possible_rank))
                possible_rank += 1

            if len(cut_points) >= n_cut_points:
                break
        #
        dict_cut_points = {idx: rank for idx, value, rank in cut_points}
        df_ent_tmp['cut-rank'] = df_ent_tmp.index.map(dict_cut_points)
        #
        # Add to list
        list_df_ent.append(df_ent_tmp)
    #
    df_ent = pd.concat(list_df_ent, axis='index')
    #
    return df_ent, df_dec


if __name__ == '__main__':

    #
    # Args
    #
    parser = argparse.ArgumentParser()
    celltypes = ['spermatocyte', 'spermatogonia', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    parser.add_argument("--celltype", default='spermatocyte', type=str, choices=celltypes, help="Cell type. Must be either 'spermatocyte' or 'enterocyte'. Defaults to spermatocyte")
    parser.add_argument("--network", default='thr', type=str, help="Network to use. Defaults to 'thr'.")
    parser.add_argument("--threshold", default=0.5, type=float, help="Threshold value. Defaults to 0.5.")
    parser.add_argument("--layer", default='DM', type=str, choices=['DM', 'MM', 'HS'], help="Network layer to compute SVD. Defaults to 'DM'.")
    parser.add_argument("--components", default=9, type=int, help="Number of singular values (components) to calculate. Defaults to 9.")
    #
    parser.add_argument("--radius_window", default=1.0, type=float, help="Window size for the radius dimension.")
    parser.add_argument("--radius_overlap", default=0.1, type=float, help="Window overlap size for the radius dimension")
    parser.add_argument("--angle_window", default=30, type=int, help="Window size for the angle dimension (in degrees)")
    parser.add_argument("--angle_overlap", default=15, type=int, help="Window overlap size for the angle dimension (in degrees)")
    #
    args = parser.parse_args()
    #
    celltype = args.celltype  # spermatocyte or enterocyte
    network = args.network
    threshold = args.threshold
    threshold_str = str(threshold).replace('.', 'p')
    layer = args.layer
    components = args.components
    #
    radius_window = args.radius_window
    radius_overlap = args.radius_overlap
    angle_window = args.angle_window
    angle_overlap = args.angle_overlap
    #
    threshold_str = str(threshold).replace('.', 'p')
    #
    #
    print('Calculating PCA entropy for {celltype:s}-{network:s}-{threshold:s}-{layer:s}'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer))
    rPCAFile = 'results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dim.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    wDiAnFile = 'results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-dian.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    wEntrFile = 'results/pca/{celltype:s}/{layer:s}/pca-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-entropy.csv.gz'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer)
    #
    df_pca = pd.read_csv(rPCAFile, index_col=0, encoding='utf-8')
    #
    df_ent, df_dian = compute_entropy(df_pca, radius_window=radius_window, radius_overlap=radius_overlap, angle_window=angle_window, angle_overlap=angle_overlap, components=9)
    #
    ensurePathExists(wDiAnFile)
    ensurePathExists(wEntrFile)
    df_ent.to_csv(wEntrFile)
    df_dian.to_csv(wDiAnFile)
