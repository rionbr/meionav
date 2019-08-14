#
# Author: Rion Brattig Correia
# Date: 24/07/2019
#
# Description: Code to compute Differential Gene Expression (DGE)
# In HS and MM we are interested in comparing:
#  - SpermatoGonia vs SpermatoCyte (interested in upregulated genes in spermatocytes)
#  - SpermatoCyte vs SpermaTid (interested in genes downregulated in spermatocytes)
#

# Necessary packages
# > install.packages("BiocManager")
# > library(BiocManager)

# Necessary Bioconductor packages
# > BiocManager::install("edgeR")
# > BiocManager::install("limma")

# Load Packages (version)
library(edgeR) # edgeR_3.26.5
library(limma) # limma_3.40.6 
library(ggplot2) # ggplot2_3.2.0

setwd('/Users/rionbr/Sites/spermnet/1-diff-gene-exp')

# Current Packages
sessionInfo()

#
# General Function to calculate DEG
#
calc_diff_gene_exp <- function(counts, group, contrasts, min.cpm, min.rowsum, file) {
  # Create DGE object
  dge <- DGEList(counts=counts, group=group)
  # Design Matrix
  formula <- formula(~0+dge$samples$group)
  design <- model.matrix(formula)
  # Rename Design Row/Cols
  rownames(design) <- colnames(dge)
  colnames(design) <- gsub("dge$samples$group", "", colnames(design), fixed=TRUE)
  # Contrasts
  contrasts <- makeContrasts(contrasts=contrasts, levels=design)
  # Remove rows with small count
  keep <- rowSums(cpm(counts) > min.cpm) >= min.rowsum
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  # trimmed mean of M- values (TMM) - affects "dge$samples$norm.factors"
  dge <- calcNormFactors(dge, method="TMMwsp")
  # Estimation
  dge <- estimateDisp(dge, design, robust=TRUE)
  ## ExatText  
  #res <- exactTest(dge)
  ### GLM Fit
  fit <- glmFit(dge, design)
  res <- glmLRT(fit, contrast=contrasts)
  ###
  # Up/Down Decision
  #status = decideTestsDGE(res, adjust.method="BH", p.value=0.01, lfc=1)
  #sumStatus <- summary(status)
  # Calculate FDR
  out <- topTags(res, adjust.method="fdr", n=Inf, sort.by="none")$table
  # Export toCSV
  write.table(out, file=file, sep=",", col.names=NA)
  # Return
  return(out)
}


####################
# [H]omo [S]apiens #
####################
cHS = read.table("data/HS_RawCounts_AllSperm.csv", header=TRUE, sep=",", row.names=1)
dHS = read.table("data/HS_Design_AllSperm.csv", header=TRUE, sep=",", row.names=1)

#
# Cyte vs Gonia
#
design <- subset(dHS, condition.high=="Cyte" | condition.high=="Gonia")
design$condition.high <- droplevels(design$condition.high)
counts <- subset(cHS, select=rownames(design))

de = calc_diff_gene_exp(
  counts=counts,
  group=design$condition.high,
  contrasts="Cyte-Gonia",
  min.cpm=1,
  min.rowsum=11,
  file="results/HS/HS-DGE_Cyte_vs_Gonia.csv")

#
# Cyte vs Tid
#
design <- subset(dHS, condition.high=="Cyte" | condition.high=="Tid")
design$condition.high <- droplevels(design$condition.high)
counts <- subset(cHS, select=rownames(design))

de = calc_diff_gene_exp(
  counts=counts,
  group=design$condition.high,
  contrasts="Tid-Cyte",
  min.cpm=1,
  min.rowsum=9,
  file="results/HS/HS-DGE_Cyte_vs_Tid.csv")

####################
# [M]us [M]usculus #
####################
cMM = read.table("data/MM_RawCounts_AllSperm.csv", header=TRUE, sep=",", row.names=1)
dMM = read.table("data/MM_Design_AllSperm.csv", header=TRUE, sep=",", row.names=1)

#
# Cyte vs Gonia
#
design <- subset(dMM, condition=="Cyte" | condition=="Gonia")
design$condition <- droplevels(design$condition)
counts <- subset(cMM, select=rownames(design))

de = calc_diff_gene_exp(
  counts=counts,
  group=design$condition,
  contrasts="Cyte-Gonia",
  min.cpm=1,
  min.rowsum=4,
  file="results/MM/MM-DGE_Cyte_vs_Gonia.csv")

#
# Cyte vs Tid
#
design <- subset(dMM, condition=="Cyte" | condition=="Tid")
design$condition <- droplevels(design$condition)
counts <- subset(cMM, select=rownames(design))

de = calc_diff_gene_exp(
  counts=counts,
  group=design$condition,
  contrasts="Tid-Cyte",
  min.cpm=1,
  min.rowsum=4,
  file="results/MM/MM-DGE_Cyte_vs_Tid.csv")

