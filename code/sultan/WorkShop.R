source("http://bioconductor.org/biocLite.R")
library(cummeRbund)
## Create a cummeRbund database from the cuffdiff output:
?readCufflinks
(cuff_data <- readCufflinks('cuffdiff_result'))
## Plot The distribution of expression for each sample:
?csDensity
# CummeRbund plot of the expression level distribution for all genes:
csDensity(genes(cuff_data))
dev.copy(pdf,file="~/Desktop/workshop1/the expression level distribution for all genes.pdf")
dev.off()
# CummeRbund plot of the expression level distribution for all isoforms:
csDensity(isoforms(cuff_data))
dev.copy(pdf,file="~/Desktop/workshop1/the expression level distribution for all isoforms.pdf")
dev.off()
# CummeRbund plot of the expression level distribution for all Transcription Start Site (TSS):
csDensity(TSS(cuff_data))
dev.copy(pdf,file="~/Desktop/workshop1/the expression level distribution for all Transcription Start Site.pdf")
dev.off()
# CummeRbund plot of the expression level distribution for all Coding Regions(CDS):
csDensity(CDS(cuff_data))
dev.copy(pdf,file="~/Desktop/workshop1/the expression level distribution for all Coding Regions(CDS).pdf")
dev.off()
## Comparing the expression of each gene in two conditions ('NIKS_Healthy', 'Infectious_EP') with a scatter plot:
csScatter(genes(cuff_data), 'NIKS_Healthy', 'Infectious_EP')
dev.copy(pdf,file="~/Desktop/workshop1/Comparing the expression of each gene in 'NIKS_Healthy', 'Infectious_EP'.pdf")
dev.off()
## Comparing the expression of each gene in two conditions ('NIKS_Healthy', 'Infectious_AA') with a scatter plot:
csScatter(genes(cuff_data), 'NIKS_Healthy', 'Infectious_AA')
dev.copy(pdf,file="~/Desktop/workshop1/Comparing the expression of each gene in 'NIKS_Healthy', 'Infectious_AA'.pdf")
dev.off()
## Comparing the expression of each gene in two conditions ('Infectious_EP', 'Infectious_AA') with a scatter plot:
csScatter(genes(cuff_data), 'Infectious_EP', 'Infectious_AA')
dev.copy(pdf,file="~/Desktop/workshop1/Comparing the expression of each gene in 'Infectious_EP', 'Infectious_AA'.pdf")
dev.off()
## Creating a volcano plot to examine differentially expressed genes:
csVolcano(genes(cuff_data), 'NIKS_Healthy', 'Infectious_EP', alpha=0.005,showSignificant=TRUE, features=T, xlimits = c(-20, 20))
dev.copy(pdf,file="~/Desktop/workshop1/a volcano plot to examine differentially expressed genes between 'NIKS_Healthy', 'Infectious_EP'.pdf")
dev.off()
## Creating a volcano plot to examine differentially expressed genes('NIKS_Healthy', 'Infectious_AA'):
csVolcano(genes(cuff_data), 'NIKS_Healthy', 'Infectious_AA',alpha=0.005,showSignificant=TRUE, features=T, xlimits = c(-10, 10))
dev.copy(pdf,file="~/Desktop/workshop1/a volcano plot to examine differentially expressed genes between 'NIKS_Healthy', 'Infectious_AA'.pdf")
dev.off()
## plotting expression levels of genes of interest with bar plots:
mygene <- getGene(cuff_data, "SLC26A2")
expressionBarplot(mygene)
expressionBarplot(isoforms(mygene))
dev.copy(pdf,file="~/Desktop/workshop1/plotting expression levels of SLC26A2 gene.pdf")
dev.off()
mygene1 <- getGene(cuff_data, "hsa-mir-6723")
expressionBarplot(mygene1)
dev.copy(pdf,file="~/Desktop/workshop1/plotting expression levels of hsa-mir-6723 gene.pdf")
dev.off()
gene_diff_data <- diffData(genes(cuff_data), 'NIKS_Healthy', 'Infectious_AA')
sin_gene_data <- subset(gene_diff_data, (significant == 'yes'))
nrow(sin_gene_data)
?cummeRbund
?csHeatmap

