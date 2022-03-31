# coding=utf-8
# Author: Rion B Correia
# Date: July 08, 2019
#
# Description: Reads GOAE results for each module and plots results
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


def plot_goea(celltype='spermatocyte', layer='DM'):

    rCSVFile = 'results/goea/goea-ortho-backbone-{layer:s}-mdlc.csv.gz'.format(layer=layer)
    df = pd.read_csv(rCSVFile)

    # Trim
    df = df.loc[(df['depth'] >= 5), :]

    # All zeros are set to the smallest computable float
    df.loc[df['p_fdr_bh'] == 0.0, 'p_fdr_bh'] = np.nextafter(0, 1)
    #
    df['1-log(p)'] = 1 - (np.log(df['p_fdr_bh']))

    print('Plotting GOEA Bars: {celltype:s} - {layer}'.format(celltype=celltype, layer=layer))

    facecolor = data[layer]['facecolor']
    specie = data[layer]['name']

    df = df.sort_values('1-log(p)', ascending=False)
    #
    dft10 = df.iloc[:10, :].sort_values('1-log(p)', ascending=True)
    sl = 75  # string slice
    dft10['name'] = dft10['name'].apply(lambda x: x[0:sl] + '..' if len(x) > sl else x)

    if len(dft10) == 0:
        print('No significant GOs.')
        return False

    # Plot
    fig, ax = plt.subplots(figsize=(4.7, 3.0))

    # P-values
    title = 'GOEA-{specie:s} ortholog backbone (mdlc mutant)'.format(specie=specie)

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
    wIMGFile = 'images/goea/img-goea-bars-{celltype:s}-{layer:s}-ortho-backbone-mdlc.pdf'.format(celltype=celltype, layer=layer)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=300, bbox_inches=None, pad_inches=0.0)
    plt.close()


def plot_wordcloud(celltype='spermatocyte', layer='DM'):

    celltype_str = celltype.title()
    rCSVFile = 'results/goea/goea-ortho-backbone-{layer:s}-mdlc.csv.gz'.format(layer=layer)
    df = pd.read_csv(rCSVFile)

    # Trim
    df = df.loc[(df['depth'] >= 5), :]

    # All zeros are set to the smallest computable float
    df.loc[df['p_fdr_bh'] == 0.0, 'p_fdr_bh'] = np.nextafter(0, 1)
    #
    df['1-log(p)'] = 1 - (np.log(df['p_fdr_bh']))
    #
    text_color = data[layer]['facecolor']
    specie = data[layer]['name']
    #
    english_stopwords = stopwords.words('english')
    print('Plotting GOEA Wordcloud: {celltype:s} - {layer}'.format(celltype=celltype, layer=layer))

    # WordCloud
    text = ' '.join(df['name'].tolist())

    if len(text) == 0:
        print('No significant GOs.')
        return False

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

    frequencies = calc_frequencies(df)
    wordcloud.generate_from_frequencies(frequencies)
    # wordcloud.generate_from_text(text)

    def color_func(word, font_size, position, orientation, random_state=None, **kwargs):
        if word in data_text_color:
            return text_color
        else:
            return 'black'
    # Recolor
    wordcloud.recolor(color_func=color_func)

    title = 'GOEA-{specie:s} ortholog backbone (mdlc mutant)'.format(specie=specie, celltype=celltype_str)
    ax.set_title(title)
    #
    wp = ax.imshow(wordcloud, interpolation='bilinear')
    #
    ax.set_xticks([])
    ax.set_yticks([])

    plt.subplots_adjust(left=0.03, right=0.97, bottom=0.17, top=0.89)
    #
    wIMGFile = 'images/goea/img-goea-wc-{celltype:s}-{layer:s}-ortho-backbone-mdlc.pdf'.format(celltype=celltype, layer=layer)
    ensurePathExists(wIMGFile)
    plt.savefig(wIMGFile, dpi=300, bbox_inches=None, pad_inches=0.0)
    plt.close()

if __name__ == '__main__':

    celltype = 'spermatocyte'  # spermatocyte or enterocyte
    #network = 'thr'  # 'thr'
    #threshold = 0.5
    #threshold_str = str(threshold).replace('.', 'p')
    layer = 'DM'

    data = {
        'HS': {
            'name': 'Human',
            'facecolor': '#d62728',
        },
        'MM': {
            'name': 'Mouse',
            'facecolor': '#7f7f7f',
        },
        'DM': {
            'name': 'Insect',
            'facecolor': '#ff7f0e',
        }
    }

    data_text_color = {'RNA', 'regulation', 'metabolic', 'process', 'protein'}

    plot_goea(celltype, layer)
    plot_wordcloud(celltype, layer)

