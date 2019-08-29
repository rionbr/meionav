# Select Core Meiotic Genes

Code in this folder is used to define which genes are shared across the three different species and therefore compose the core meiotic genes.


Instructions:
1. Run `map_de_genes_to_string.py` to match DE genes (from previous step folder) to STRING-DB identifiers.
2. Run `compute_core_genes.py` to retrieve the genes with homolog overlaps (using EggNOG).
3. Run `generate_testing_genes.py` to generate the list of genes that can be tested in each different species.


Other files:
- `plot_diff_gene_exp.py` plots the results from previous step (`../1-diff-gene-exp/`) while including information about the core genes.