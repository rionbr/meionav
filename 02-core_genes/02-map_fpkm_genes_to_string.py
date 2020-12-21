# coding=utf-8
# Author: Rion B Correia
# Date: Jul 05, 2019
#
# Description: Maps FPKM Spermatogonia genes to String-DB
#
import tarfile
import pandas as pd
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
from utils import open_undefined_last_column_files
from pybiomart import Dataset


if __name__ == '__main__':

    data_replicates = {
        'spermatogonia': {
            'HS': 11,
            'MM': 3,
            'DM': 2,
        },
        'spermatocyte': {
            'HS': 12,
            'MM': 3,
            'DM': 2,
        },
        'spermatid': {
            'HS': 6,
            'MM': 3,
            'DM': 2,
        },
        'enterocyte': {
            'HS': 3,
            'MM': 3,
            'DM': 2,
        },
        'neuron': {
            'HS': 4,
            'MM': 3,
            'DM': 6,
        },
        'muscle': {
            'HS': 6,
            'MM': 6,
            'DM': 1,
        }
    }

    # Load String Files
    print("Load String File")
    data_string = {
        'HS': open_undefined_last_column_files('../data/StringDB/9606/9606.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"]),
        'MM': open_undefined_last_column_files('../data/StringDB/10090/10090.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"]),
        'DM': open_undefined_last_column_files('../data/StringDB/7227/7227.protein.aliases.v11.0.txt.gz', skiprows=1, n_fixed_cols=2, names=["id_string", "alias", "source"]),
    }

    # Load BioMart Data
    print("Load BioMart data")
    data_biomart = {
        'HS': Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org').query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID'),
        'MM': Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org').query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype', 'description']).set_index('Gene stable ID'),
        'DM': Dataset(name='dmelanogaster_gene_ensembl', host='http://www.ensembl.org').query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID'),
    }
    celltype_DM_translate = {
        'spermatogonia': 'apical',
        'spermatocyte': 'middle',
        'spermatid': 'basal'
    }

    for celltype in ['spermatogonia', 'spermatocyte', 'spermatid', 'enterocyte', 'neuron', 'muscle']:
        print('Celltype: {celltype:s}'.format(celltype=celltype))

        for specie in ['HS', 'MM', 'DM']:
            print('Specie: {specie:s}'.format(specie=specie))

            # DM has a different celltype name (Apical, Middle & Basal)
            if (specie == 'DM') and (celltype in ['spermatogonia', 'spermatocyte', 'spermatid']):
                celltypet = celltype_DM_translate.get(celltype)
            else:
                celltypet = celltype

            # Load FPKM tar.gz file
            tfile = '../01-diff-gene-exp/data/FPKM/{specie:s}/{specie:s}_{celltype:s}.tar.gz'.format(specie=specie, celltype=celltypet.title())
            tar = tarfile.open(tfile, "r:gz")

            # Replicates
            number_of_replicates = data_replicates[celltype][specie]

            # Extract and Load DataFrames
            ldf = []
            for replicate in range(1, number_of_replicates + 1):
                file = '{specie:s}_{celltype:s}_rep{replicate:d}.txt'.format(specie=specie, celltype=celltype.title(), replicate=replicate)
                print(file)
                tdf = pd.read_csv(tar.extractfile(file), sep='\t', index_col='Gene ID', usecols=['Gene ID', 'FPKM', 'TPM'])
                ldf.append(tdf)
            df = pd.concat(ldf, axis='index', sort=True)

            # Mean FPKM
            dfg = df.groupby(df.index).mean()
            dfg.index.name = 'id_gene'

            # Remove numbers from gene id
            dfg.index = dfg.index.map(lambda x: x.split('.')[0])

            # Parse String Data - Note some genes have multiple id_string, others have no match
            df_SA = data_string[specie]
            df_SA = df_SA.loc[df_SA['alias'].isin(dfg.index.to_list()), ["alias", "id_string"]].rename(columns={"alias": "id_gene"})
            df_SAg = df_SA.groupby('id_gene').agg({'id_string': lambda x: x if len(x) == 1 else list(x)})
            dfg['id_string'] = df_SAg['id_string']

            # Map BioMart
            df_BM = data_biomart[specie]
            dfg['gene'] = df_BM['Gene name']
            dfg['biotype'] = df_BM['Gene type']

            # Index Rows/Cols
            maskcols = ['id_string', 'gene', 'FPKM', 'TPM', 'biotype']
            dfg = dfg.loc[:, maskcols]
            print(dfg.head())

            # To CSV
            wCSVfile = 'results/FPKM/{specie:s}/{specie:s}-FPKM-{celltype:s}.csv.gz'.format(specie=specie, celltype=celltype)
            dfg.to_csv(wCSVfile)

    print('Done.')
