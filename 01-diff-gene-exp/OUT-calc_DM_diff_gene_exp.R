#
# Author: Rion Brattig Correia
# Date: 24/07/2019
#
# Description: Code to compute Differential Gene Expression (DGE) from MicroArray data.
# In DM we are interested in comparing:
#  - Mid Testis va Apical (interested in upregulated genes in Mid)
#  - Mid vs Distal Testis (interested in genes downregulated in Mid)
#

# Necessary packages
# > install.packages("BiocManager")
# > library(BiocManager)

# Necessary Bioconductor packages
# > BiocManager::install("affy")
# > BiocManager::install("limma")

# Load Packages (version)
library(affy) # affy_1.62.0
library(limma) # limma_3.40.6
library(ggplot2) # ggplot2_3.2.0

setwd('/Users/rionbr/Sites/spermnet/1-diff-gene-exp')

# Current Packages
sessionInfo()

#
# General Function to calculate DEG
# In comparison to HS and MM: reads == counts;
#
calc_diff_gene_exp <- function(reads, group, contrasts, file) {
  # Design Matrix
  formula <- formula(~0+group)
  design <- model.matrix(formula)
  # Rename Design Row/Cols
  rownames(design) <- colnames(reads)
  colnames(design) <- gsub("group", "", colnames(design), fixed=TRUE)
  # Contrasts
  contrasts <- makeContrasts(contrasts=contrasts, levels=design)
  # Estimation
  fit <- lmFit(reads, design)
  fit <- contrasts.fit(fit, contrasts)
  res <- eBayes(fit)
  # Calculate FDR
  out <- topTable(res, coef=1, number=nrow(res), adjust.method='BH', sort.by="none")
  # Export toCSV
  write.table(out, file=file, sep=",", col.names=NA)
  # Return
  return(out)
}

###############################
# [D]rosophila [M]elanogaster #
###############################
dDM = read.table("data/DM/DM_Design_Raw.csv", header=TRUE, sep=",", row.names=1)

celfile.path = "data/DM/GSE18502_RAW"
# Load cel.gz Files
data <- ReadAffy(celfile.path=celfile.path,  sampleNames=rownames(dDM), verbose=TRUE)
# Phenotypic Data
pData(data) <- dDM
# Normalize Data
eset = rma(data, normalize=TRUE, bgversion=2, verbose=TRUE)
#eset <- mas5(data)
# Calls (P/M/A)
calls <- mas5calls(data)
# To DataFrame
cDM = data.frame(exprs(eset))


#
# Mid vs Apical
#
design <- subset(dDM, condition=="Mid" | condition=="Apical")
design$condition <- droplevels(design$condition)
reads <- subset(cDM, select=rownames(design))

de = calc_diff_gene_exp(
  reads=reads,
  group=design$condition,
  contrasts="Mid-Apical",
  file="results/DGE/DM/DM-DGE_Mid_vs_Apical-(MicroArray).csv")

#
# Mid vs Distal
#
design <- subset(dDM, condition=="Mid" | condition=="Distal")
design$condition <- droplevels(design$condition)
reads <- subset(cDM, select=rownames(design))

de = calc_diff_gene_exp(
  reads=reads,
  group=design$condition,
  contrasts='Mid-Distal',
  file="results/DGE/DM/DM-DGE_Mid_vs_Distal-(MicroArray).csv")
