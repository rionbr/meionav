# coding=utf-8
# Author: Rion B Correia
# Date: May 19, 2020
#
# Description: Calculates the number of expressed genes (protein-coding & non-protein coding) for all cell types across different TPM cut-offs
#
import pandas as pd


if __name__ == '__main__':

    species = ['HS', 'MM', 'DM']
    celltypes = ['spermatogonia', 'spermatocyte', 'spermatid', 'enterocyte', 'neuron', 'muscle']
    thresholds = [0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]

    r = []
    for specie in species:

        print("Calculating for species: {specie:s}".format(specie=specie))

        for celltype in celltypes:

            print("Calculating for celltype: {celltype:s}".format(celltype=celltype))

            rFPKMfile = 'results/{specie:s}-FPKM-{celltype:s}.csv.gz'.format(specie=specie, celltype=celltype)
            df = pd.read_csv(rFPKMfile)
            print(df.head())

            for threshold in thresholds:

                dft = df.loc[df['TPM'] >= threshold, :]
                n_genes = len(dft)
                biotype_counts = dft['biotype'].value_counts()
                n_pc_genes = biotype_counts.loc['protein_coding']
                n_not_pc_genes = biotype_counts.loc[biotype_counts.index != 'protein_coding'].sum()
                #
                r.append((specie, celltype, threshold, n_genes, n_pc_genes, n_not_pc_genes))

    dfR = pd.DataFrame(r, columns=['species', 'celltype', 'threshold', 'genes', 'pc-genes', 'not-pc-genes'])

    #
    # Export
    #
    dfR.to_csv('results/comparison-FPKM-thresholds.csv.gz')
