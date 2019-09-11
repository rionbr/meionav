# Select Core Meiotic Genes

Code in this folder is used to define which genes are shared across the three different species and therefore compose the core meiotic genes. In practice, there are multiple ways one could define what are the 'core' genes, thus we have different pipelines for each assessment.


Instructions:
1. Run `map_de_genes_to_string.py` to match DE genes (from previous step folder) to STRING-DB identifiers.
2. Run `compute_core_genes.py` to select from EggNOG only those meta-genes with associated ids among our species.
2. Run `compute_core_genes-<pipeline>.py` to then subselect meta-genes with homolog overlaps. Results are then saved withint each `/results/<pipeline-name>` subfolder.


Note: Currently there are two possible pipelines: `conserved` and `pooling`, both which are computed for different values of FDR.

Other files:
- `plot_diff_gene_exp.py` plots the results from previous step (`../1-diff-gene-exp/`) while including information about the core genes.