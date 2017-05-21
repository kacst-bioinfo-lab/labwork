source("http://bioconductor.org/biocLite.R")
library("GenomicFeatures")
library("Rsamtools")
fls <- list.files("tophat_all/", pattern="bam$", full.names =T)
bamfiles <- BamFileList(fls)
(sampleTable <- read.csv("sample_table.csv",row.names = 1))
## Check the chromosome names using seqinfo function
seqinfo(bamfiles[1])
##Defining Gene Model
##Reading the gene model from Ensembl GTF file using makeTxDbFromGFF from the GenomicFeatures package.
## I want to make a list of exons grouped by gene for counting read/fragments.
gtffile <- ("Homo_sapiens.GRCh37.75.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character()))
##producing a GRangesList of all exons grouped by gene
(ebg <- exonsBy(txdb, by="gene"))
## Read Counting Step
library("GenomicAlignments")
se <- summarizeOverlaps(features = ebg, reads = bamfiles, mode = "Union", singleEnd= F, ignore.strand = T, fragments= T)
str(metadata(rowRanges(se)))
colData(se)
(colData(se) <- DataFrame(sampleTable))
## The DESeqDataSet, Sampleinformation, and designFormula.
##checking the millions of fragments that uniquely aligned to the genes.
round(colSums(assay(se))/1e6, 1)
biocLite("DESeq2")
library("DESeq2")
se$condition <- relevel(se$condition, "Healthy_NIKS")
dds <- DESeqDataSet(se, design = ~ condition)
countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)
ddMat <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
## Pre-filtering the Data:
nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds)
## the rlog transformation:
rld <- rlog(dds, blind = F)
head(assay(rld), 3)
## Showing effect of transformation
par( mfrow = c(1,2))
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized = T)[,1:2] + 1), pch=16 , cex=0.3)
plot(assay(rld)[,1:2], pch=16, cex=0.3)
## Sample distances:
sampleDists <- dist(t(assay(rld)))
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.copy(pdf,file="~/Desktop/workshop1/Heatmap of sample to sample distance using rlog-transformed values.pdf")
dev.off()
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd ) 
rownames(samplePoisDistMatrix) <- paste( rld$condition) 
colnames(samplePoisDistMatrix) <- NULL 
pheatmap(samplePoisDistMatrix,
                        clustering_distance_rows=poisd$dd,
                        clustering_distance_cols=poisd$dd,
                        col=colors)
## PCA Plot 2D 
plotPCA(rld, intgroup = c("condition"))
dev.copy(pdf,file="~/Desktop/workshop1/PCA plot using the rlog-transformed values.pdf")
dev.off()
(data <- plotPCA(rld, intgroup = c( "condition", "SampleName"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
ggplot(data, aes(PC1, PC2, color=SampleName, shape=condition)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.copy(pdf,file="~/Desktop/workshop1/PCA plot2 using the rlog-transformed values.pdf")
dev.off()
dds <- DESeq(dds)
(res <- results(dds))
mcols(res, use.names=TRUE)
summary(res)
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("condition"))
plotMA(res, ylim=c(-4,4))
