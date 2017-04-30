#####################################################################
############### Title : RNA-Seq WorkFlow analysis Using Bioconductor#
############### Author: Sultan Alharbi ##############################
############### Date: 27/04/2017 ####################################
#####################################################################
## Load package
library("airway")
data("airway")
se <- airway
se$dex <- relevel(se$dex, "untrt")
## I will quickly check the millions of fragments that uniquely aligned to genes. the second argument of round tells how many decimal points to keep) 
round(colSums(assay(se)) / 1e6, 1)
colData(se)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)
## Starting from count matrices
countdata <- assay(se)
(head(countdata, 3))
coldata <- colData(se)
## countdata: a table with the fragment counts.
## coldata: a table with information about the samples.
## To now construct the DESeqDataSet object from the matrix of counts and the sample information table:
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ cell + dex))
## Pre-filtering the dataset
## Checking the number of rows (genes).
(nrow(dds))
## I will remove rows of DESeqDataSet that have no counts or only a single count across all samples:
((dds <- dds[rowSums(counts(dds)) > 1, ]))
(nrow(dds))
## normalization the data
(rld <- rlog(dds, blind = F))
head(assay(rld), 3)
## blind=F means that differences between cell lines and treatment should not add to the variance-mean profile of the experiment
## the effect of the transformation, I will plot the first sample against the sceond usning the log2 function
(par(mfrow = c(1, 2)))
(dds <- estimateSizeFactors(dds))
## Scatterplot using the log2
(plot(log2(counts(dds, normalized=T)[,1:2] +1), pch=16, cex=0.3)) 
## Scatterplot using the rlog
(plot(assay(rld)[,1:2], pch=16, cex=0.3))
## Sample distances.
## overall simalirity between Samples
## using dist function to calculate the euclidean distance between samples.
(sampleDists <- dist(t(assay(rld))))
### visualizing the distances in heatmap.
library("pheatmap")
library("RColorBrewer")
(sampleDistMatrix <- as.matrix(sampleDists))
(rownames(sampleDistMatrix) <- paste(rld$dex, rld$cell, sep="-"))
(colnames(sampleDistMatrix) <- NULL)
(colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255))
(pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col= colors))
library("PoiClaClu")
(poisd <- PoissonDistance(t(counts(dds))))
(samplePoisDisMatrix <- as.matrix(poisd$dd))
(rownames(samplePoisDisMatrix) <- paste(rld$dex, rld$cell, sep="-"))
(colnames(samplePoisDisMatrix) <- NULL)
(pheatmap (samplePoisDisMatrix, clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col= colors))
## Princioal Components Analysis (PCA) plot.
(plotPCA(rld, intgroup= c("dex","cell")))
