### BACKUP ###

#
# [H]omo [S]apiens
#

# Load Data
cHS = read.table("data/HS_RawCounts_AllSperm.csv", header=TRUE, sep=",", row.names=1)
dHS = read.table("data/HS_Design_AllSperm.csv", header=TRUE, sep=",", row.names=1)

# Only keep the following libraries:
# > Spermatogonia_Ad_1, Spermatogonia_Ad_2,Spermatogonia_Ad_4, Spermatogonia_Ad_6
# > Spermatogonia_Ap_1, Spermatogonia_Ap_2, Spermatogonia_Ap_3, Spermatogonia_Ap_4, Spermatogonia_Ap_5
# > Spermatocytes_Early_1, Spermatocytes_Early_3, Spermatocytes_Early_4, Spermatocytes_Early_5, Spermatocytes_Early_6
# > Spermatocytes_Late_2
keep = c(1,2,4,6,7,8,9,10,11,12,14,15,16,17,19,24,25,26,27,28,29)
cHSt = cHS[, keep]
dHSt = dHS[keep ,]


gtUp = read.table('goldstandard/HS_Upregulated_CytesVSGonia_edgeR.csv', header=TRUE, sep=',', row.names=1, quote="\"")
gtUp = rownames(gtUp)
gtDown = read.table('goldstandard/HS_Downregulated_CytesVSTids_edgeR.csv', header=TRUE, sep=',', row.names=1, quote="\"")
gtDown = rownames(gtDown)

#
# HS - Spermatocytes vs Spermatogonia
#
# Check Libraries
d <- subset(dHS, condition.high=="Gonia" | condition.high=="Cyte")
d$condition.high <- droplevels(d$condition.high)
c <- subset(cHS, select=rownames(d))
dge <- DGEList(counts=c, group=d$condition.high)
pseudoCounts <- log2(dge$counts+1)
plotMDS(pseudoCounts, top=500, dim.plot=c(1,2), cex=1, gene.selection='pairwise', col=c(rep("darkred",6),rep("red",5),rep("darkblue",6),rep("blue",6)),)
#
# Remove outlier libraries
remove = c("Spermatocytes_Early_2","Spermatocytes_Early_3")
cHSt = cHS[, -which(names(cHS) %in% remove)]
dHSt = dHS[-which(names(cHS) %in% remove) ,]
#
# Recompute DGE
d <- subset(dHSt, condition.high=="Cyte" | condition.high=="Gonia")
d$condition.high <- droplevels(d$condition.high)
c <- subset(cHSt, select=rownames(d))
# Create DGE object
dge <- DGEList(counts=c, group=d$condition.high)
# Design Matrix
formula <- formula(~0+dge$samples$group)
design <- model.matrix(formula)
# Rename Design Row/Cols
rownames(design) <- colnames(dge)
colnames(design) <- gsub("dge$samples$group", "", colnames(design), fixed=TRUE)
# Contrasts
contrasts <- makeContrasts(contrasts="Cyte-Gonia", levels=design)
# Remove rows with small count
keep <- rowSums(cpm(c) > 1) >= 11
dge <- dge[keep, , keep.lib.sizes=FALSE]
# trimmed mean of M- values (TMM) - affects "dge$samples$norm.factors"
dge <- calcNormFactors(dge, method="TMM")
# Estimation
dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmFit(dge, design)
res <- glmLRT(fit, contrast=contrasts)
# Up/Down Decision
status = decideTestsDGE(res, adjust.method="BH", p.value=0.01, lfc=1)
sumStatus <- summary(status)
# Export toCSV
out <- topTags(res, adjust.method="BH", n=Inf, sort.by="none")$table
write.table(out, file="result/HS-DGE_Cyte_vs_Gonia.csv", sep=",", col.names=NA)
write.table(out[status==1,], file="result/HS-DGE_UpCyte_vs_Gonia.csv", sep=",", col.names=NA)
# Plot
limma::plotMD(res, status=status, main=paste("MD Plot:"), hl.col=alpha(c("firebrick", "blue"), 0.4), values=c(1, -1), xlab="Average Expression", ylab="logFC")
abline(h=0, col="grey", lty=2)

tUp = rownames(subset(status, status==1))

length(gtUp)
length(tUp)

length(intersect(gtUp, tUp))
length(intersect(gtUp, tUp)) / length(union(gtUp, tUp))
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
length(outersect(gtUp,tUp))

plotBCV(dge, main="BCV Plot")
pseudoCounts <- log2(dge$counts+1)
plotMDS(pseudoCounts, col=c(rep("blue",11), rep("red",10)),)


#
# HS - Spermatids vs Spermatocytes
#
d <- subset(dHS, condition.high=="Cyte" | condition.high=="Tid")
d$condition.high <- droplevels(d$condition.high)
c <- subset(cHS, select=rownames(d))
# Create DGE object
dge <- DGEList(counts=c, group=d$condition.high)
# Design Matrix
formula <- formula(~0+dge$samples$group)
design <- model.matrix(formula)
# Rename Design Row/Cols
rownames(design) <- colnames(dge)
colnames(design) <- gsub("dge$samples$group", "", colnames(design), fixed=TRUE)
# Contrasts
contrasts <- makeContrasts(contrasts="Tid-Cyte", levels=design)
# Remove rows with small count
keep <- rowSums(cpm(c) > 1) >= 9
dge <- dge[keep, , keep.lib.sizes=FALSE]
# trimmed mean of M- values (TMM) - affects "dge$samples$norm.factors"
dge <- calcNormFactors(dge, method="TMM")
# Estimation
dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmFit(dge, design)
res <- glmLRT(fit, contrast=contrasts)
# Up/Down Decision
status = decideTestsDGE(res, adjust.method="BH", p.value=0.01, lfc=1)
sumStatus <- summary(status)
# Export toCSV
out <- topTags(res, adjust.method="BH", n=Inf, sort.by="none")$table
write.table(out, file="result/HS-DGE_Cyte_vs_Tid.csv", sep=",", col.names=NA)
write.table(out[status==-1,], file="result/HS-DGE_DownCyte_vs_Tid.csv", sep=",", col.names=NA)
# Plot
limma::plotMD(res, status=status, main=paste("MD Plot:"), hl.col=alpha(c("firebrick", "blue"), 0.4), values=c(1, -1), xlab="Average Expression", ylab="logFC")
abline(h=0, col="grey", lty=2)

tDown = rownames(subset(status, status==-1))

length(gtDown)
length(tDown)

length(intersect(gtDown, tDown))
length(intersect(gtDown, tDown)) / length(union(gtDown, tDown))
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
length(outersect(gtDown,tDown))


plotBCV(dge, main="BCV Plot")
pseudoCounts <- log2(dge$counts+1)
plotMDS(pseudoCounts)

#
# [M]us [M]usculus
#
# Load Data
cMM = read.table("data/MM_RawCounts_AllSperm.csv", header=TRUE, sep=",", row.names=1)
dMM = read.table("data/MM_Design_AllSperm.csv", header=TRUE, sep=",", row.names=1)

gtDown = read.table('goldstandard/MM_MeioticGenes_DownRegulated.csv', header=TRUE, sep=',', row.names=1, quote="\"")
gtDown = rownames(gtDown)
gtUp = read.table('goldstandard/MM_MeioticGenes_UpRegulated.csv', header=TRUE, sep=',', row.names=1, quote="\"")
gtUp = rownames(gtUp)

#
# MM - Spermatocytes vs Spermatogonia
#
d <- subset(dMM, condition=="Cyte" | condition=="Gonia")
d$condition <- droplevels(d$condition)
c <- subset(cMM, select=rownames(d))
# Create DGE object
dge <- DGEList(counts=c, group=d$condition)
# Design Matrix
formula <- formula(~0+dge$samples$group)
design <- model.matrix(formula)
# Rename Design Row/Cols
rownames(design) <- colnames(dge)
colnames(design) <- gsub("dge$samples$group", "", colnames(design), fixed=TRUE)
# Contrasts
contrasts <- makeContrasts(contrasts="Cyte-Gonia", levels=design)
# Remove rows with small count
keep <- rowSums(cpm(c) > 1) >= 4
dge <- dge[keep, , keep.lib.sizes=FALSE]
# trimmed mean of M- values (TMM) - affects "dge$samples$norm.factors"
dge <- calcNormFactors(dge, method="TMM")
# Estimation
dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmFit(dge, design)
res <- glmLRT(fit, contrast=contrasts)
# Up/Down Decision
status = decideTestsDGE(res, adjust.method="BH", p.value=0.01, lfc=1)
sumStatus <- summary(status)
# Export toCSV
out <- topTags(res, adjust.method="BH", n=Inf, sort.by="none")$table
write.table(out, file="result/MM-DGE_Cyte_vs_Gonia.csv", sep=",", col.names=NA)
write.table(out[status==1,], file="result/MM-DGE_UpCyte_vs_Gonia.csv", sep=",", col.names=NA)

tUp = rownames(subset(status, status==1))
length(gtUp)
length(tUp)
length(intersect(gtUp, tUp))
length(intersect(gtUp, tUp)) / length(union(gtUp, tUp))
#
# MM - Spermatids vs Spermatocytes
#
d <- subset(dMM, condition=="Tid" | condition=="Cyte")
d$condition <- droplevels(d$condition)
c <- subset(cMM, select=rownames(d))
# Create DGE object
dge <- DGEList(counts=c, group=d$condition)
# Design Matrix
formula <- formula(~0+dge$samples$group)
design <- model.matrix(formula)
# Rename Design Row/Cols
rownames(design) <- colnames(dge)
colnames(design) <- gsub("dge$samples$group", "", colnames(design), fixed=TRUE)
# Contrasts
contrasts <- makeContrasts(contrasts="Tid-Cyte", levels=design)
# Remove rows with small count
keep <- rowSums(cpm(c) > 1) >= 4
dge <- dge[keep, , keep.lib.sizes=FALSE]
# trimmed mean of M- values (TMM) - affects "dge$samples$norm.factors"
dge <- calcNormFactors(dge, method="TMM")
# Estimation
dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmFit(dge, design)
res <- glmLRT(fit, contrast=contrasts)
# Up/Down Decision
status = decideTestsDGE(res, adjust.method="BH", p.value=0.01, lfc=1)
sumStatus <- summary(status)
# Export toCSV
out <- topTags(res, adjust.method="BH", n=Inf, sort.by="none")$table
write.table(out, file="result/MM-DGE_Cyte_vs_Tid.csv", sep=",", col.names=NA)
write.table(out[status==-1,], file="result/MM-DGE_DownCyte_vs_Tid.csv", sep=",", col.names=NA)

tDown = rownames(subset(status, status==-1))
length(gtDown)
length(tDown)
length(intersect(gtDown, tDown))
length(intersect(gtDown, tDown)) / length(union(gtDown, tDown))
# Note: one gene difference between old and new pipeline: ENSMUSG00000025758

