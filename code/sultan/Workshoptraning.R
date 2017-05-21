source("http://bioconductor.org/biocLite.R")
biocLite(c("pathview", "gage", "gageData", "GenomicAlignments"))
##Reading the gene model from Ensembl GTF file using makeTxDbFromGFF from the GenomicFeatures package.
## I want to make a list of exons grouped by gene for counting read/fragments.
library("GenomicFeatures")
gtffile <- ("Homo_sapiens.GRCh37.75.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character()))
exByGn <- exonsBy(txdb, "gene")
fls <- list.files("/Users/sultanalharbi/SULTAN_Dropbox/SULTAN_BIN/tophat_all/", pattern="bam$", full.names =T)
bamfiles <- BamFileList(fls)
flag <- scanBamFlag(isSecondaryAlignment=FALSE, isProperPair=TRUE)
param <- ScanBamParam(flag=flag)
#to run multiple core option: library(parallel); options("mc.cores"=8)
library(parallel)
options("mc.cores"=4)
gnCnt <- summarizeOverlaps(exByGn, bamfiles, mode = "Union", ignore.strand = T, singleEnd=F, param=param)
cnts=assay(gnCnt)
dim(cnts)
## Genes with 0 counts across all samples are removed.
sel.rn=rowSums(cnts) !=0
cnts=cnts[sel.rn,]
dim(cnts)
## I will divide the read counts by the total number of mapped reads for each sample as to normalize over library size and sequence depth.
libsizes=colSums(cnts)
size.factor=libsizes/exp(mean(log(libsizes)))
cnts.norm=t(t(cnts)/size.factor)
range(cnts.norm)
## I will add a small yet appropriate positive count (+8) to all genes before doing log2 transformation as to avoid -inf
cnts.norm=log2(cnts.norm+8)
range(cnts.norm)
##optional MA plot.
?par
?plot
?abline
pdf("cnts.maplots.pdf", width=8, height=10)
op=par(lwd=2, cex.axis=1.5, cex.lab=1.5, mfrow=c(2,1))
plot((cnts.norm[,7]+cnts.norm[,5])/2, (cnts.norm[,7]-cnts.norm[,5]), main="(a) Normal vs Normal", xlab="mean", ylab="change",ylim=c(-5,5), xlim=c(0,20), lwd=1)
abline(h=0, lwd=2, col="red", lty="dashed")
plot((cnts.norm[,1]+cnts.norm[,5])/2, (cnts.norm[,1]-cnts.norm[,5]), main="(b) Infectious vs Infectious", xlab="mean", ylab="change", ylim=c(-5,5), xlim=c(0,20), lwd=1)
abline(h=0, lwd=2, col="red", lty="dashed")
dev.off()
library(gage)
?gage
head(kegg.gs)
ref.idx=1:2
samp.idx=3:7
data("kegg.gs")
cnts.kegg.pre <- gagePrep(cnts.norm, ref = ref.idx, samp = samp.idx, same.dir = TRUE, compare = "as.group", rank.test = FALSE, use.fold = TRUE, weights = NULL, full.table =FALSE)
cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref= ref.idx, samp = samp.idx, compare = "as.group", rank.test = FALSE, use.fold = TRUE, saaPrep = gagePrep)
cnts.d= cnts.norm[, samp.idx]-rowMeans(cnts.norm[, ref.idx])
sel <- cnts.kegg.p$greater[, "q.val"] < 0.1 &
  +         !is.na(cnts.kegg.p$greater[,"q.val"])
path.ids <- rownames(cnts.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8)
library(pathview)
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = cnts.d, pathway.id = pid, species = "hsa"))
sel.l <- cnts.kegg.p$less[, "q.val"] < 0.1 &
  +            !is.na(cnts.kegg.p$less[,"q.val"])
path.ids.l <- rownames(cnts.kegg.p$less)[sel.l]
path.ids.l2 <- substr(path.ids.l, 1, 8)
pv.out.list.l <- sapply(path.ids.l2[1:3], function(pid) pathview(gene.data = cnts.d, pathway.id = pid, species = "hsa"))
