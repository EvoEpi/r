library(gplots)
library(tidyr)

setwd("/PATH/TO/YOUR/WORKING/DIRECTORY")

data.1<-read.table("botrytis_cinerea.apothecium_disk.abun.out",
                   header=TRUE,
                   sep="\t",
                   check.names=FALSE,
                   comment.char="")
data.1$Tissue<-"Apothecium disk"

data.2<-read.table("botrytis_cinerea.apothecium_primordia.abun.out",
                   header=TRUE, sep="\t",
                   check.names=FALSE,
                   comment.char="")
data.2$Tissue<-"Apothecium primordia"

data.3<-read.table("botrytis_cinerea.apothecium_stipes.abun.out",
                   header=TRUE,
                   sep="\t",
                   check.names=FALSE,
                   comment.char="")
data.3$Tissue<-"Apothecium stipes"

data.total<-rbind(data.1, data.2, data.3)

data.sub<-select(data.total, "Gene ID", Tissue, FPKM)

data.spread<-spread(data.sub, key=Tissue, value=FPKM)

#replace numbered row names with gene names
rownames(data.spread)<-data.spread$"Gene ID"

#remove genes name column, which is the first column
data.spread[,1]<-NULL

#convert table to a matrix...a matrix must have the same data type (numeric,
#character, etc.) and the same length. A data frame is more general than a
#matrix, in that different columns can have different modes (numeric, character,
#factor, etc.)
matrix<-as.matrix(data.spread)

#log10 transform all values in matrix
log.matrix<-log10(matrix+1)

#a function to calculate variance per row...i.e., per gene variance across tissues
rowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

#pass your matrix through the variance function
var.genes<-rowVar(log.matrix)

#sort genes based on variance and get the top 500 most variable genes
top.genes<-names(sort(var.genes, decreasing=TRUE))[1:500]

#add expression values to most variable genes
var.fpkm<-log.matrix[top.genes,]

#create a color palette for plotting
my.palette<-colorRampPalette(c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                "#e0f3f8", "#fee090", "#fdae61", "#f46d43",
                                "#d73027", "#a50026"))(n=50)

#saving the plot to a pdf is a little different...you must first declare a file
#to save to...
pdf("botrytis_cinerea.tissue.abun.100_most_var.heatmap.pdf", height=5, width=7)
#make your heatmap and send it to a variable h.map
h.map<-heatmap.2(var.fpkm, scale="column", col=my.palette, trace="none",
                tracecol="#000000", cexRow=0.1, cexCol=1, margins=c(8,2),
                srtCol=45)
#close the R environment to let R know you are done adding things to the pdf
dev.off()
