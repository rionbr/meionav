from data_enterocyte import *
from data_spermatocyte import *

#
# Placeholder
#
default_layer_placeholder = {
    'HS': None,
    'MM': None,
    'DM': None}

#
# Enterocyte
#
enterocyte_modules_pca = {
    'modules': {
        'HS': enterocyte_modules_pca_hs,
        'MM': enterocyte_modules_pca_mm,
        'DM': enterocyte_modules_pca_dm
    },
    'files': dict(default_layer_placeholder),
    'dfs': dict(default_layer_placeholder)}

data_enterocyte = {
    'graphs': dict(default_layer_placeholder),
    'modules-pca': enterocyte_modules_pca}

#
# Spermatocyte
#
spermatocyte_modules_pca = {
    'modules': {
        'HS': spermatocyte_modules_pca_hs,
        'MM': spermatocyte_modules_pca_mm,
        'DM': spermatocyte_modules_pca_dm
    },
    'files': dict(default_layer_placeholder),
    'dfs': dict(default_layer_placeholder)}

data_spermatocyte = {
    'graphs': dict(default_layer_placeholder),
    'modules-pca': spermatocyte_modules_pca}

#
# Cells
#
data_cells = {
    'spermatocyte': data_spermatocyte,
    'enterocyte': data_enterocyte
}
