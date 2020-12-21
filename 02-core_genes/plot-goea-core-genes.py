# coding=utf-8
# Author: Rion B Correia
# Date: Nov 17, 2020
#
# Description: Reads GOAE results for core genes and plots results
#
#
import math
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import ensurePathExists
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['mathtext.fontset'] = 'cm'
import matplotlib.pyplot as plt
#
from wordcloud import WordCloud
from nltk.corpus import stopwords


def plot_goea(df, celltype='spermatocyte', species='HS', facecolor='red', ns='BP'):

    df = df.copy()
    # Select
    df = df.loc[(df['NS'] == ns), :]
    # Trim
    df = df.loc[(df['depth'] >= 5), :]

    # All zeros are set to the smallest computable float
    df.loc[df['p_fdr_bh'] == 0.0, 'p_fdr_bh'] = np.nextafter(0, 1)
    #
    df['1-log(p)'] = 1 - (np.log(df['p_fdr_bh']))

    print('Plotting GOEA Bars: {celltype:s} {species} {ns:s}'.format(celltype=celltype, species=species, ns=ns))
    species_str = dict_species[species]
    ns_str = dict_ns[ns]

    df = df.sort_values('1-log(p)', ascending=False)
    #
    dft10 = df.iloc[:10, :].sort_values('1-log(p)', ascending=True)
    sl = 75 # string slice
    dft10['name'] = dft10['name'].apply(lambda x: x[0:sl] + '..' if len(x) > sl else x)

    if len(dft10) == 0:
        print('No significant GOs.')
        return None

    # Plot
    fig, ax = plt.subplots(figsize=(4.7, 3.0))

    # P-values
    title = 'GO enrichment - {species:s} {ns:s}'.format(species=species_str, ns=ns_str)

    ind = np.arange(0, len(dft10), 1)
    bp = ax.barh(ind, 1 - np.log(dft10['p_fdr_bh']), height=0.8, facecolor=facecolor, zorder=4)
    ax.set_title(title, fontsize='large')

    minx, maxx = ax.get_xlim()
    for bar, name in zip(bp.patches, dft10['name'].tolist()):
        bx = bar.get_x()
        by = bar.get_y()
        bh = bar.get_height()
        # bw = bar.get_width()
        tx = bx + (0.01 * maxx)
        ty = (by + (bh / 2))
        ax.text(x=tx, y=ty, s=name, ha='left', va='center', fontsize='x-small', zorder=5)
    #
    ax.axvline(x=(1 - math.log(0.01)), color='#666666', ls='dotted')
    ax.axvline(x=(1 - math.log(0.05)), color='#c7c7c7', ls='dashed')
    ax.set_yticks(ind)
    ax.set_yticklabels(dft10['GO'])
    ax.set_xlabel(r'$1 - $log($p$-value)')
    ax.set_ylim(-0.7, (10 - 0.3))
    ax.grid(axis='x', zorder=1)

    plt.subplots_adjust(left=0.21, right=0.97, bottom=0.17, top=0.89)
    #plt.tight_layout()
    #
    wIMGFile = 'images/goea-bars/img-goea-bars-{celltype:s}-{species:s}-core-genes-{ns:s}.pdf'.format(celltype=celltype, species=species, ns=ns)
    print(wIMGFile)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=300, bbox_inches=None, pad_inches=0.0)
    plt.close()


def plot_wordcloud(df, celltype='spermatocyte', species='DM', facecolor='red', ns='BP'):

    celltype_str = celltype.title()

    # Trim
    df = df.loc[(df['depth'] >= 5), :]

    # All zeros are set to the smallest computable float
    df.loc[df['p_fdr_bh'] == 0.0, 'p_fdr_bh'] = np.nextafter(0, 1)
    #
    df['1-log(p)'] = 1 - (np.log(df['p_fdr_bh']))

    species_str = dict_specie[species]
    #
    english_stopwords = stopwords.words('english')
    print('Plotting GOEA Wordcloud: {celltype:s} - {species}'.format(celltype=celltype, species=species))

    # WordCloud
    dft = df.loc[(df['module-id'] == mid), :]
    text = ' '.join(dft['name'].tolist())

    if len(text) == 0:
        print('No significant GOs.')
        return None

    text = text.replace('-', ' ')
    #
    fig, ax = plt.subplots(figsize=(4.0, 3.0))

    def color_func(*args, **kwargs):
        return (0, 0, 0)

    wordcloud = WordCloud(background_color='white', max_font_size=45, width=400, height=300, stopwords=english_stopwords, relative_scaling='auto', colormap='tab10', color_func=color_func, collocation_threshold=20)

    def calc_frequencies(dfA):
        r = []
        for i, dfAt in dfA.iterrows():
            name = dfAt['name']
            pvalue = dfAt['1-log(p)']
            name = name.replace('-', ' ').replace(',', '').replace('.', '').replace("'", '')
            for word in name.split(' '):
                if word not in english_stopwords:
                    r.append((i, word, pvalue))

        dfr = pd.DataFrame(r, columns=['id', 'name', 'pvalue']).set_index('id')
        dfr['name'] = dfr['name'].replace('proteasomal', 'proteasome')
        #
        dfrg = dfr.groupby('name').agg({'pvalue': ['count', 'sum']})
        dfrg.columns = dfrg.columns.droplevel()
        dfrg['frequency'] = dfrg['count'].rank(method='min') * dfrg['sum'].rank(method='min')
        dfrg.sort_values('frequency', ascending=False, inplace=True)

        return dfrg.reset_index().set_index('name')['frequency'].to_dict()

    frequencies = calc_frequencies(dft)
    wordcloud.generate_from_frequencies(frequencies)
    # wordcloud.generate_from_text(text)

    def color_func(word, font_size, position, orientation, random_state=None, **kwargs):
        if word in data_text_color[mid]:
            return text_color
        else:
            return 'black'
    # Recolor
    wordcloud.recolor(color_func=color_func)

    title = 'GOEA-{specie:s} {celltype:s} M{mid:d}-{mname:s}'.format(specie=specie, celltype=celltype_str, mid=mid, mname=dict_replace[mname])
    ax.set_title(title)
    #
    wp = ax.imshow(wordcloud, interpolation='bilinear')
    #
    ax.set_xticks([])
    ax.set_yticks([])

    plt.subplots_adjust(left=0.03, right=0.97, bottom=0.17, top=0.89)
    #
    wIMGFile = 'images/goea-wordcloud/{celltype:s}/{layer:s}/img-goea-wc-{celltype:s}-{network:s}-{threshold:s}-{layer:s}-mod-{mid:d}.pdf'.format(celltype=celltype, network=network, threshold=threshold_str, layer=layer, mid=mid)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=300, bbox_inches=None, pad_inches=0.0)
    plt.close()


if __name__ == '__main__':

    celltype = 'spermatocyte'  # spermatocyte or enterocyte
    species = 'HS'

    dict_species = {'HS': 'Human', 'MM': 'Mouse', 'DM': 'Insect'}
    dict_ns = {'BP': 'biological processes', 'MF': 'molecular functions'}

    """
    data_text_color = {
        1: ['ubiquitination', 'ubiquitin'],
        2: ['splicing'],
        3: ['translation', 'translational', 'cotranslational'],
        4: ['rRNA'],
        5: ['vesicle', 'transport'],
        6: ['respiration', 'respiratory', 'electron'],
        7: ['cell', 'cycle'],
        8: ['DNA', 'repair'],
        9: ['mitochondrial', 'translation', 'translational'],
        10: ['cell', 'cycle'],
        11: ['metabolic'],
        12: ['histidine', 'peptidyl', 'dephosphorylation'],
    }
    """
    for species in ['HS', 'MM', 'DM']:

        rCSVFile = 'results/goea/goea-{celltype:s}-{species:s}-core-genes.csv.gz'.format(celltype=celltype, species=species)
        df = pd.read_csv(rCSVFile)
        print(df['NS'].value_counts())

        plot_goea(df, celltype, species, facecolor='#2ca02c', ns='BP')
        plot_goea(df, celltype, species, facecolor='#98df8a', ns='MF')
        #
        #plot_wordcloud(df, celltype, species, facecolor='#d62728', ns='BP')
        #plot_wordcloud(df, celltype, species, facecolor='#1f77b4', ns='MF')
