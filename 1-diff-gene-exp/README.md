# Calculates Differential Gene Expression (DGE)

Code in this folder is used to define which genes are Differentially Expressed.


Instructions:
1. Run `calc_diff_gene_exp.R` to compute gene expression levels (logFC, FDR and p-value) for all genes in HS and MM.
2. Run `select_dge_genes.py` to select only those genes that are DGE and that we are interested in. Note two different FDR levels are used (FDR<=0.05 and 0.01).


Other files:
- `plot_DM_diff_gene_exp.py` and `plot_HS_MM_diff_gene_exp.py` plots the results from (1).
- `calc_plot_HS_PCA.py` and `calc_plot_MM_PCA.py` calculate Principal Component Analysis (PCA) and plot the results.

## R Session Info
```R
> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.2.0 edgeR_3.26.5  limma_3.40.6 

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1       rstudioapi_0.10  magrittr_1.5     splines_3.6.1    tidyselect_0.2.5
 [6] munsell_0.5.0    statmod_1.4.32   colorspace_1.4-1 lattice_0.20-38  R6_2.4.0        
[11] rlang_0.4.0      dplyr_0.8.3      tools_3.6.1      grid_3.6.1       gtable_0.3.0    
[16] withr_2.1.2      lazyeval_0.2.2   assertthat_0.2.1 tibble_2.1.3     crayon_1.3.4    
[21] purrr_0.3.2      glue_1.3.1       compiler_3.6.1   pillar_1.4.2     scales_1.0.0    
[26] locfit_1.5-9.1   pkgconfig_2.0.2 
```

---

## Original RNAseq data


Mus Musculus:
- [GSE43717](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43717)

> Soumillon et al. *Cell Rep*. 2013 Jun 27;3(6):2179-90. doi: [10.1016/j.celrep.2013.05.031](https://www.sciencedirect.com/science/article/pii/S2211124713002489)



Homo Sapiens:
- [SRP069329](https://www.ncbi.nlm.nih.gov/sra/?term=SRP069329)

> Jan, et al. *Development* 2017 144: 3659-3673; doi: [10.1242/dev.152413](https://dev.biologists.org/content/144/20/3659)