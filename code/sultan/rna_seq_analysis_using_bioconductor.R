#####################################################################
############### Title : RNA-Seq WorkFlow analysis Using Bioconductor#
############### Author: Sultan Alharbi ##############################
############### Date: 27/04/2017 ####################################
#####################################################################
## Load package
source("https://bioconductor.org/biocLite.R")
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
##Showing the effect of the tranformation
## Scatterplot using the log2
(plot(log2(counts(dds, normalized=T)[,1:2] +1), pch=16, cex=0.3)) 
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/scatterplot_using_log2.pdf")
dev.off()
## Scatterplot using the rlog
(plot(assay(rld)[,1:2], pch=16, cex=0.3))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/scatterplot_using_rlog.pdf")
dev.off()
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
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/visualizing_distances_between_samples_heatmap.pdf")
dev.off()
library("PoiClaClu")
(poisd <- PoissonDistance(t(counts(dds))))
(samplePoisDisMatrix <- as.matrix(poisd$dd))
(rownames(samplePoisDisMatrix) <- paste(rld$dex, rld$cell, sep="-"))
(colnames(samplePoisDisMatrix) <- NULL)
(pheatmap (samplePoisDisMatrix, clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col= colors))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/visualizing_distances_between_samples_heatmap1.pdf")
dev.off()
## Princioal Components Analysis (PCA) plot.Sample to Sample distance (2D plane)
(plotPCA(rld, intgroup= c("dex","cell")))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/PCA0.pdf")
dev.off()
(data <- plotPCA(rld, intgroup=c("dex", "cell"), returnData=T))
(percentVar <- round(100 * attr(data, "percentVar")))
library("ggplot2")
(ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) + xlab(paste0("PC1:", percentVar[1], "% variance" )) + ylab(paste0 ("PC2:", percentVar[2], "% variance")))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/PCA1.pdf")
dev.off()
## Runing the differential expression pipeline
(dds <- DESeq(dds))
(res <- results(dds))
(mcols(res, use.names = T))
(summary(res))
(res.05 <- results(dds, alpha = .05))
(table(res.05$padj < .05))
(topGene <- rownames(res)[which.min(res$padj)])
(plotCounts(dds, gene=topGene, intgroup = c("dex")))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/normalized_counts_for_singlegene1.pdf")
dev.off()
(data <- plotCounts(dds, gene=topGene, intgroup=c("dex","cell"), returnData=TRUE))
(ggplot(data, aes(x=dex, y=count, color=cell)) + scale_y_log10() + geom_point(position=position_jitter(width=.1,height=0), size=3))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/normalized_counts_for_singlegene2.pdf")
dev.off()
(ggplot(data, aes(x=dex, y=count,  fill=dex)) + scale_y_log10() + geom_dotplot(binaxis="y", stackdir="center"))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/normalized_counts_for_singlegene3.pdf")
dev.off()
(ggplot(data, aes(x=dex, y=count, color=cell, group=cell)) + scale_y_log10() + geom_point(size=3) + geom_line())
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/normalized_counts_for_singlegene4.pdf")
dev.off()
(plotMA(res, ylim=c(-5,5)))
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/MA_plot.pdf")
dev.off()
(resLFC1 <- results(dds, lfcThreshold = 1))
## Labelling individual points on the MA-plot
(plotMA(resLFC1, ylim=c(-5,5)))
(topGene <- rownames(resLFC1)[which.min(resLFC1$padj)])
with(resLFC1[topGene, ], {
    points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2) 
    text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")})
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/labelling_significant_genes.pdf")
dev.off()
## Histogram of the p values.
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20, col="grey50", border = "white")
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/histogram_of_the_pvalues.pdf")
dev.off()
## Gene Clustering
## selecting the 20 genes that have the highest variance across samples.
library("genefilter")
(topVarGenes <- head(order(rowVars(assay(rld)), decreasing = T), 20))
mat <- assay(rld)[topGene, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)
## Independent filtering_For weakly expressed genes_ratio of small p values by normalized count.
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2))) 
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE)) 
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
dev.copy(pdf,file="/Users/sultanalharbi/OneDrive/labwork/data/sultan/theratio_samll_pvalue.pdf")
dev.off()
## annotating and exporting results.
library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
res$entrez <- mapIds(org.Hs.eg.db, keys = row.names(res), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
resOrdered <- res[order(res$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF,  file="results.csv")
library("ReportingTools")
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <-  finish(htmlRep)
browseURL(url)
(resGR <- results(dds, lfcThreshold=1, format="GRanges"))
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")
library("Gviz")
window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)
sig <- factor(ifelse(resGRsub$padj < .1 & !is.na(resGRsub$padj),"sig","notsig"))
options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig) 
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
                type="h", name="log2 fold change", strand="+") 
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")
