# coding=utf-8
# Author: Rion B Correia
# Date: May 19, 2020
#
# Description: Calculates the number of expressed genes (protein-coding & non-protein coding) for all cell types across different TPM cut-offs
#
import numpy as np
import pandas as pd
from pybiomart import Dataset


if __name__ == '__main__':

    species = ['HS', 'MM', 'DM']
    celltypes = ['spermatogonia', 'spermatocyte', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    thresholds = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]
    #thresholds = np.arange(0.1, 9, 0.5)

    #
    datamart_names = {'HS': 'hsapiens_gene_ensembl', 'MM': 'mmusculus_gene_ensembl', 'DM': 'dmelanogaster_gene_ensembl'}

    r = []
    for specie in species:

        print("Calculating for species: {specie:s}".format(specie=specie))

        print("Querying Datamart")
        datamart_name = datamart_names[specie]
        ds = Dataset(name=datamart_name, host='http://www.ensembl.org')
        dfQ = ds.query(attributes=['ensembl_gene_id', 'external_gene_name', 'gene_biotype']).set_index('Gene stable ID')
        #
        n_genome = len(dfQ)
        dfQpc = dfQ.loc[(dfQ['Gene type'] == 'protein_coding'), :]
        n_genome_pc = len(dfQpc)
        DFQnpc = dfQ.loc[(dfQ['Gene type'] != 'protein_coding'), :]
        n_genome_non_pc = len(DFQnpc)
        print('done.')

        for celltype in celltypes:

            print("Calculating for celltype: {celltype:s}".format(celltype=celltype))

            rFPKMfile = '../../02-core_genes/results/FPKM/{specie:s}/{specie:s}-FPKM-{celltype:s}.csv.gz'.format(specie=specie, celltype=celltype)
            df = pd.read_csv(rFPKMfile)

            for threshold in thresholds:

                dft = df.loc[df['TPM'] >= threshold, :]
                n_genes = len(dft)
                biotype_counts = dft['biotype'].value_counts()
                n_genes_pc = biotype_counts.loc['protein_coding']
                n_genes_non_pc = biotype_counts.loc[biotype_counts.index != 'protein_coding'].sum()
                #
                r.append((specie, celltype, threshold, n_genome, n_genome_pc, n_genome_non_pc, n_genes, n_genes_pc, n_genes_non_pc))
    #
    dfR = pd.DataFrame(r, columns=['species', 'celltype', 'threshold', 'genome', 'genome-pc', 'genome-non-pc', 'genes', 'genes-pc', 'genes-non-pc'])

    #
    # Export
    #
    dfR.to_csv('results/species-celltypes-threshold-stats.csv.gz')
