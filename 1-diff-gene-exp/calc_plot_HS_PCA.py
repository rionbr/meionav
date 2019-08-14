# coding=utf-8
# Author: Rion B Correia
# Date: Jul 24, 2019
#
# Description: Plots HS and MM genes that are differentially expressed, calculated by 1_calc_diff_gene_exp.R
#
# Instructions:
#   For each species, two different comparisons are necessary:
#   Spermatocytes vs Spermatogonia (interested in genes upregulated in spermatocytes)
#   Spermatids vs Spermatocytes (interested in genes downregulated in spermatids)

#   Genes are considered upregulated if: log2FC >1 +  Pval < 0.01 + FDR≤0.05
#   Genes are considered downregulated if: log2FC <-1 +  Pval < 0.01 + FDR≤0.05
#
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances


def labelScatterPlot(labels, ax, xs, ys, n, stdn, printLabel):
    df = pd.DataFrame({'xs':xs,'ys':ys}, index=labels)
    stdx , stdy = df['xs'].std() , df['ys'].std()
    df = df.loc[ (df['xs'] >= stdx*stdn) | (df['xs'] <= -(stdx*stdn)) | (df['ys'] >= stdy*stdn) | (df['ys'] <= -(stdy*stdn)) ]

    df = df.sort_values('xs', ascending=False)
    textobjs = list()
    for label, data in df.iterrows():
        textobj = ax.text(data['xs'], data['ys'], label, fontsize=8, alpha=1, zorder=10, ha='center', va='center')
        textobjs.append(textobj)
    return textobjs


def annotatePoint(label, ax, x, y, xpad, ypad):
    obj = ax.annotate(
        label, 
        xy = (x, y), xytext = (x, y),
        xycoords='data', textcoords='data', ha='center', va='center',
        bbox = dict(boxstyle = 'round,pad=0.2', fc = 'yellow', alpha = 0.3),
        #arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0.5', alpha = 0.2),
        fontsize=9, alpha=0.65
        )
    return obj


if __name__ == '__main__':

    dfC = pd.read_csv('data/HS_RawCounts_AllSperm.csv', index_col=0, nrows=None).T
    dfD = pd.read_csv('data/HS_Design_AllSperm.csv', index_col=0, nrows=None)

    #dfC = pd.read_csv('data/MM_RawCounts_AllSperm.csv', index_col=0, nrows=None).T
    #dfD = pd.read_csv('data/MM_Design_AllSperm.csv', index_col=0, nrows=None)

    # Loc
    dfD = dfD.loc[ (dfD['condition-high'].isin(['Spermatogonia','Spermatocytes'])), : ]
    dfC = dfC.loc[ (dfD.index), : ]

    remove = ['Spermatocytes_Early_6','Spermatogonia_Ap_5','Spermatogonia_Ad_6','Spermatocytes_Late_6','Spermatocytes_Early_3','Spermatocytes_Early_2','Spermatocytes_Late_5','Spermatogonia_Ap_4','Spermatogonia_Ap_3']
    remove = ['Spermatocytes_Early_6']
    dfD = dfD.loc[ ~(dfD.index.isin(remove)), : ]
    dfC = dfC.loc[ ~(dfC.index.isin(remove)), : ]

    # PseudoCount
    pc = np.log2(dfC.values + 1)

    print('> Calculating Distances')
    D = pairwise_distances(pc, metric='euclidean')
    #D_clipped = np.clip(D, a_min=-1, a_max=1200)

    print('> Plot Distances')
    fig, ax1 = plt.subplots(figsize=(7,7),  nrows=1, ncols=1)

    hm = ax1.imshow(D, cmap='jet_r', interpolation='nearest')
    
    # Separating Lines
    #for xy in [5.5,10.5,16.5,22.5]:
    #    ax1.axhline(xy, lw=0.5, color='black')
    #    ax1.axvline(xy, lw=0.5, color='black')

    ax1.set_title('HS Library Euclidean Distance')
    ax1.set_yticks(range(dfD.shape[0]))
    ax1.set_xticks(range(dfD.shape[0]))
    ax1.set_yticklabels(dfD.index.tolist())
    ax1.set_xticklabels(dfD.index.tolist(), ha='center', va='top', rotation=90)
    plt.colorbar(hm, label='distance')

    plt.gca().invert_yaxis()    

    plt.tight_layout()
    #plt.subplots_adjust(left=0.07, right=0.98, bottom=0.06, top=0.95, wspace=0.32, hspace=0.35)
    plt.savefig('images/img-HS-dist.pdf', dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()

    #asd

    print('> Normalizing')
    X = dfC.values
    X = StandardScaler().fit_transform(X)


    print('> PCA')
    pca = PCA(n_components=None)
    res = pca.fit_transform(X)

    print('> to DataFrame')
    dfPCA = pd.DataFrame(res[:, 0:9], columns=['1c', '2c', '3c', '4c', '5c', '6c', '7c', '8c', '9c'], index=dfC.index)
    dfPCA['label'] = dfC.index
    dfPCA['condition-high'] = dfD['condition-high']
    dfPCA['condition-low'] = dfD['condition-low']
    dfPCA['color-high'] = dfPCA['condition-high'].map({
        'Spermatocytes':'#d62728',
        'Spermatogonia':'#1f77b4'
        })
    dfPCA['color-low'] = dfPCA['condition-low'].map({
        'Spermatogonia_Ad':'#d62728',
        'Spermatogonia_Ap':'#ff9896',
        'Spermatocytes_Early':'#1f77b4',
        'Spermatocytes_Late':'#aec7e8',
        'Spermatids':'green',
        })

    s = pd.Series(pca.explained_variance_ratio_, index=range(1, res.shape[1] + 1))

    print('> Plotting')

    # Files
    #
    title = r'HS PCA'    

    fig, ((ax1,ax2,ax3),(ax4,ax5,ax6), (ax7,ax8,ax9)) = plt.subplots(figsize=(12,10),  nrows=3, ncols=3)

    fig.suptitle(title)
    # plt.rc('font', size=11)
    # plt.rc('figure', titlesize=12)
    # plt.rc('legend', fontsize=12)
    # plt.rc('legend', scatterpoints=1)

    #
    # Plot EigenVals
    #
    s_cumsum = s.cumsum()
    n_eigen_95 = s_cumsum[s_cumsum<0.95].shape[0]

    n = 9
    ind = np.arange(n)
    height = s.iloc[:n].values
    width = 0.60
    xticklabels = ind+1

    cmap = mpl.cm.get_cmap('hsv_r')
    norm = mpl.colors.Normalize(vmin=0,vmax=n)
    s_colors = list(map(cmap, np.linspace(0,1,n, endpoint=False)))

    print(ind)
    print(height)
    print(width)
    print(s_colors)
    ax1.bar(ind, height, width, color=s_colors, zorder=9, edgecolor='black', align='center')
    ax1.set_xticks(ind)
    ax1.set_xticklabels(xticklabels)

    ax1.set_title('Explained variance ratio')
    ax1.annotate('95% with {:d}\nsingular vectors'.format(n_eigen_95), xy=(0.97, 0.97), xycoords="axes fraction", ha='right', va='top')
    ax1.set_xlabel('Components')
    ax1.set_ylabel('%')

    ax1.grid()
    ax1.set_xlim(-.5,n-.5)

    for dim, ax in zip( range(1,10), [ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9] ):
        print('> Dim: {:d}'.format(dim,dim+1))
        col = str(dim)+'c'
        x = str(dim)+'c'
        y = str(dim+1)+'c'
        xs = dfPCA[x].tolist()
        ys = dfPCA[y].tolist()
        #pca_colors = dfPCA['color'].tolist()

        #if pca_colors == None:
        #   pca_colors = ys
        #   norm = mpl.colors.Normalize( vmin=min(ys), vmax=max(ys) )
        #   ax.scatter(xs,ys, c=pca_colors, cmap=cmapR2G2G, norm=norm, marker='o', edgecolor='black', lw=0.5, s=25)
        #else:
        ax.scatter(xs, ys, c=dfPCA['color-low'], marker='o', edgecolor='black', lw=0.5, s=30, zorder=5, rasterized=True)
        
        # Draw a X at the center
        #ax.plot(0,0, color='black', marker='x', ms=16)
        
        # Draw lines at the center
        ax.axhline(y=0, c='gray', lw=0.75, ls='-', zorder=2)
        ax.axvline(x=0, c='gray', lw=0.75, ls='-', zorder=2)

        ax.set_title('Components {} and {}'.format(dim, dim+1) )
        ax.set_xlabel('Component %d' % (dim))
        ax.set_ylabel('Component %d' % (dim+1))
        
        
        ax.grid()
        ax.axis('equal')
        ax.locator_params(axis='both', tight=True, nbins=6)
        
        #ax.set_aspect('equal')
        #xmin, xmax = ax.get_xlim()
        #ax.set_xticks(np.round(np.linspace(xmin, xmax, 6), 2))
        #ymin, ymax = ax.get_ylim()
        #ax.set_yticks(np.round(np.linspace(ymin, ymax, 6), 2))
        
        labels = labelScatterPlot(dfPCA['label'].values, ax, xs, ys, n=2, stdn=2, printLabel=True)
        """
        adjust_text(
            labels, x=xs, y=ys, ax=ax,
            lim=1000,
            force_text=(.1, .4), force_points=(.1, .4), force_objects=(1, 1),
            expand_text=(1.4, 1.8), expand_points=(1.4, 1.8), expand_objects=(1,1), expand_align=(1.1,1.9),
            only_move={'points':'xy','text':'xy','objects':'xy'},
            text_from_points=True,
            ha='center',va='center', autoalign=False,
            arrowprops=dict(arrowstyle="->", shrinkB=5, color='gray', lw=0.5, connectionstyle='angle3'),
        )
        """
        

    # Save
    #plt.tight_layout()
    plt.subplots_adjust(left=0.07, right=0.98, bottom=0.06, top=0.95, wspace=0.32, hspace=0.35)
    plt.savefig('images/img-HS-pca.pdf', dpi=150, bbox_inches=None, pad_inches=0.0)
    plt.close()
