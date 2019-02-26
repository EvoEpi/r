#install libraries if necessary, e.g., install.packages("ggplot2")
#load libraries
library(ggplot2)
library(reshape2)
library(dplyr)

setwd("/PATH/TO/YOUR/WORKING/DIRECTORY")

#read in gene expression adundance files from hisat2
data.1<-read.table("botrytis_cinerea.apothecium_disk.abun.out",
                   header=TRUE,
                   sep="\t",
                   check.names=FALSE,
                   comment.char="")
#add a column specifying tissue library was made from
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

#join data frames (datasets) vertically
data.total<-rbind(data.1, data.2, data.3)

#select columns you want
data.sub<-select(data.total, "Gene ID", Tissue, FPKM)

#histogram
#log10(FPKM+1)=add 1 to FPKM values and then log10 transform
ggplot(data.sub, aes(x=log10(FPKM+1), y=..density.., fill=Tissue)) +
  #each plus (+) at the end of a line represents a new command for ggplot2
  #bins=Number of bins
  #position=Position adjustment, either as a string, or the result of a call to a position adjustment function
  geom_histogram(bins=10, position="dodge") +
  #manually change colour of bars
  scale_fill_manual(values=c("#E69F00", "#0072B2", "#CC79A7")) +
  #change x-axis label
  xlab("Gene expression (log10(FPKM+1))") +
  #change y-axis label
  ylab("Density") +
  #theme=Control all non-data display
  theme_classic()

#density
ggplot(data.sub, aes(x=log10(FPKM+1), fill=Tissue)) +
  #alpha=transparency
  geom_density(alpha=0.25) +
  #facet_wrap wraps a 1d sequence of panels into 2d
  facet_wrap(Tissue~.) +
  xlab("Gene expression (log10(FPKM+1))") +
  ylab("Density") +
  theme_classic()

#box plot
ggplot(data.sub, aes(x=Tissue, y=log10(FPKM+1))) +
  geom_boxplot(notch=TRUE) +
  xlab("Tissue") +
  ylab("Gene expression (log10(FPKM+1))") +
  theme_classic()

#swap FPKM=0 with FPKM=NA
data.sub$FPKM[data.sub$FPKM==0]<-NA

#density with NAs removed
ggplot(data.sub, aes(x=log10(FPKM+1), fill=Tissue)) +
  #na.rm removes NA
  geom_density(alpha=0.25, na.rm=TRUE) +
  xlab("Gene expression (log10(FPKM+1))") +
  ylab("Density") +
  theme_classic()
