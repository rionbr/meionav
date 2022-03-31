import os
import gzip
from io import StringIO
import pandas as pd


def get_network_layer(G, layer=''):
    return G.subgraph([n for n, d in G.nodes(data=True) if (d.get('layer') == layer)]).copy()


def open_undefined_last_column_files(filepath, skiprows=0, n_fixed_cols=None, *args, **kwargs):
    """ Some StringDB files need manual parsing to be loaded as a pandas DataFrame."""
    with gzip.open(filepath, 'rt') as f:
        ios = u''
        # Skip header
        for i in range(skiprows):
            _ = f.readline()
        # Loop file content
        for i, line in enumerate(f, start=0):
            sline = line.split('\t')
            ios += u'\t'.join(sline[:n_fixed_cols]) + u'\t' + sline[-1]

        return pd.read_csv(StringIO(ios), sep='\t', encoding='utf-8', *args, **kwargs)


def ensurePathExists(path):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
