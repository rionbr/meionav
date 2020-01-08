# Select Core Meiotic Genes

Code in this folder is used to define which genes are shared across the three different species and therefore compose the core meiotic genes. In practice, there are multiple ways one could define what are the 'core' genes. We opted for two different pipelines: `mammals' and `core'. The first only considers HS and MM genes, the second includes genes also present in DM.


Instructions:
1. Run `map_de_genes_to_string.py` to match DE genes (from previous step folder) to STRING-DB identifiers.
2. Run `calc_meta_genes.py` to select from EggNOG only those meta-genes with associated ids among our species.
2. Run `calc_pipeline-<pipeline>.py` to then subselect meta-genes with homolog overlaps. Results are then saved withint each `/results/<pipeline-name>` subfolder.


Other files:
- `plot_diff_gene_exp.py` plots the results from previous step (`../1-diff-gene-exp/`) while including information about the core genes.