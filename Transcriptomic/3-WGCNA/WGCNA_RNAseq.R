###############################################################################
# Script: WGCNA_RNAseq.R
# Description:
#   Perform WGCNA on normalized RNA-seq expression data and relate modules 
#   to environmental traits.
###############################################################################

# Libraries (install once if needed)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("WGCNA", "edgeR", "cluster", "DESeq2"))
# install.packages(c("flashClust","data.table"))

library(WGCNA)
library(flashClust)
library(data.table)
library(DESeq2)

# Set working directory
setwd("/path/to/WGCNA_directory")

# -------------------------
# STEP 1: Load and preprocess data
# -------------------------

# Load counts
counts <- read.delim("counts.txt")
rownames(counts) <- make.unique(as.character(counts[,1]))
counts <- counts[,-1]
counts <- round(counts)

# Filter lowly expressed genes (mean count >= 0.5)
cat("Genes before filtering:", nrow(counts), "\n")
keep_expr <- rowMeans(counts) >= 0.5
counts <- counts[keep_expr, ]
cat("Genes after filtering:", nrow(counts), "\n")

# DESeq2 object with dummy condition
fake_conditions <- data.frame(row.names = colnames(counts), condition = factor(rep("X", ncol(counts))))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = fake_conditions, design = ~1)

# VST normalization
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Filter low variance genes
variances <- apply(vsd_mat, 1, var)
fixed_threshold <- 0.3
keep_variance <- variances >= fixed_threshold
vsd_mat <- vsd_mat[keep_variance, ]
cat("Genes after variance filtering:", nrow(vsd_mat), "\n")

# Z-score normalization
vsd_mat_zscore <- t(scale(t(vsd_mat)))
write.table(vsd_mat_zscore, file = "counts_norm_zscore.txt", sep = "\t", quote = FALSE)

# Load expression data for WGCNA
datExpr <- read.table("counts_norm_zscore.txt", sep = "\t", header = TRUE)
datExpr <- as.data.frame(t(datExpr))
rownames(datExpr) <- gsub("^X", "", rownames(datExpr))

# Check for sample/gene outliers
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

# Load traits
datTraits <- read.delim2("parasitisme.txt")
rownames(datTraits) <- datTraits$Condition
datTraits$Condition <- NULL
stopifnot(all(rownames(datTraits) == rownames(datExpr)))

# Cluster samples and plot dendrogram
sampleTree <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

# -------------------------
# STEP 2: Network construction and module detection
# -------------------------

allowWGCNAThreads()
powers <- c(1:10, seq(12, 40, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed hybrid")

# Plot scale-free topology and mean connectivity
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.90, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, col="red")

# Use chosen soft threshold
softPower <- 8
adjacency <- adjacency(datExpr, power = softPower, type = "signed hybrid")
TOMadj <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOMadj

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main="Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# Module detection
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels=FALSE, hang=0.03, main="Gene dendrogram and module colors")

# Merge similar modules
dynamic_MEList <- moduleEigengenes(datExpr, colors=dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes
dynamic_MEDiss <- 1 - cor(dynamic_MEs)
merge_thresh <- 0.2
merge_dynamic <- mergeCloseModules(datExpr, dynamicColors, cutHeight = merge_thresh, verbose = 3)

moduleColors <- merge_dynamic$colors
mergedMEs <- merge_dynamic$newMEs

gene2module <- data.frame(Gene = names(datExpr), Module = moduleColors)
write.csv(gene2module, "Gene2Module.csv", row.names = FALSE)
write.csv(mergedMEs, "MEs.csv")

# -------------------------
# STEP 3: Relate modules to traits
# -------------------------
nSamples <- nrow(datExpr)
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
names(MEs) <- substring(names(MEs), 3)

moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep="")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text=0.6, zlim=c(-1,1), main="Module-trait Relationships")

# -------------------------
# STEP 4: Gene significance & module membership
# -------------------------
trait <- as.data.frame(datTraits$Pythium)
names(trait) <- "trait"
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use="p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
geneTraitSignificance <- as.data.frame(cor(datExpr, trait, use="p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste("GS.", names(trait), sep="")
names(GSPvalue) <- paste("p.GS.", names(trait), sep="")

geneInfo <- data.frame(Gene = names(datExpr), Module = moduleColors, GS = geneTraitSignificance[,1], MM = NA)
for(i in 1:nrow(geneInfo)) {
  mod <- geneInfo$Module[i]
  if(mod %in% names(geneModuleMembership)) geneInfo$MM[i] <- geneModuleMembership[i, mod]
}
write.csv(geneInfo, "GS-MM_trait.csv", row.names = FALSE)

# Extract candidate genes with high GS & MM
FilterGenes <- abs(geneTraitSignificance) > 0.7 & abs(geneModuleMembership$salmon) > 0.7
interesting_genes <- names(datExpr)[FilterGenes]
write.csv(interesting_genes, "Genes_trait.csv")
