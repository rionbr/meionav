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
enterocyte_modules_svd = {
    'modules': {
        'HS': enterocyte_modules_svd_hs,
        'MM': enterocyte_modules_svd_mm,
        'DM': enterocyte_modules_svd_dm
    },
    'files': dict(default_layer_placeholder),
    'dfs': dict(default_layer_placeholder)}

data_enterocyte = {
    'graphs': dict(default_layer_placeholder),
    'modules-svd': enterocyte_modules_svd}

#
# Spermatocyte
#
spermatocyte_modules_svd = {
    'modules': {
        'HS': spermatocyte_modules_svd_hs,
        'MM': spermatocyte_modules_svd_mm,
        'DM': spermatocyte_modules_svd_dm
    },
    'files': dict(default_layer_placeholder),
    'dfs': dict(default_layer_placeholder)}

data_spermatocyte = {
    'graphs': dict(default_layer_placeholder),
    'modules-svd': spermatocyte_modules_svd}

#
# Cells
#
data_cells = {
    'spermatocyte': data_spermatocyte,
    'enterocyte': data_enterocyte
}
