#!/usr/bin/env Rscript

require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(ChIPseeker)
require(clusterProfiler)
require(ggplot2)
require(UpSetR)
require(org.Dm.eg.db)

#######
#######
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
chip_dmel<- read.table("counts_table.bed")
colnames(chip_dmel)<- c("chrom","start","end","name","pileup", "strand", "fold", "-log10(pvalue)","-log10(qvalue)","chip_dme","input_dme")
head(chip_dmel)

plot1<-ggplot(chip_dmel) + geom_point(aes(x=log2(chip_dme +1 ), y=log2(input_dme +1))) +
  geom_abline(intercept = 0, slope = 1, col="red")

print(plot1)
dev.copy(pdf, file="Figures/counts.pdf")
dev.off()

gr2<-GRanges(seqnames = chip_dmel$chrom, ranges = IRanges(chip_dmel$start,chip_dmel$end))
head(gr2)

peak<-readPeakFile("short_peaks.narrowPeak")
peak

plot2<-covplot(peak, weightCol="X35")

print(plot2)

dev.copy(pdf, file="Figures/covplot.pdf")
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=100, downstream=100)
tagMatrix <- getTagMatrix(peak, windows=promoter)

plot3<-tagHeatmap(tagMatrix, xlim=c(-100, 100), color="red")

print(plot3)

dev.copy(pdf, file="Figures/tagHeatmap.pdf")
dev.off()

plot4<-peakHeatmap(peak, TxDb=txdb, upstream=100, downstream=100, color="red")

print(plot4)

dev.copy(pdf, file="Figures/peakHeatmap.pdf")
dev.off()

plot5<-plotAvgProf(tagMatrix, xlim=c(-100, 100), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

print(plot5)

dev.copy(pdf, file="Figures/plotAvgProf.pdf")
dev.off()


plot6<-plotAvgProf2(peak, TxDb=txdb, upstream=100, downstream=100, 
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
print(plot6)

dev.copy(pdf, file="Figures/plotAvgProf2.pdf")
dev.off()

plot7<-plotAvgProf(tagMatrix, xlim=c(-100, 100), conf = 0.95, resample = 1000)

print(plot7)

dev.copy(pdf, file="Figures/tagAvgProf.pdf")
dev.off()
#####
peakAnno <- annotatePeak(gr2,
                         tssRegion=c(-1000,1000), 
                         TxDb = txdb,
                         annoDb="org.Dm.eg.db")
as.GRanges(peakAnno)
#######################
plot8<-plotAnnoPie(peakAnno)
print(plot8)

dev.copy(pdf, file="Figures/AnnoPie.pdf")
dev.off()


plot9<-vennpie(peakAnno)
print(plot9)
dev.copy(pdf, file="Figures/VennPie.pdf")
dev.off()


plot10<-plotAnnoBar(peakAnno)
print(plot10)

dev.copy(pdf, file="Figures/AnnotBar.pdf")
dev.off()


plot11<-upsetplot(peakAnno)
print(plot11)

dev.copy(pdf, file="Figures/upsetplot.pdf")
dev.off()


q()
