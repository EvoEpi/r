library(ggplot2); library(dplyr)

setwd("/YOUR/WORKING/DIRECTORY")

#infiles are bedtools intersect files of a three column "makewindows" bed file and four column feature bed file
d.gene<-read.table("bt-intersect.CorFlo.w100000_s50000-highQgm.gene.tsv",
                   sep="\t",
                   header=FALSE,
                   check.names=FALSE,
                   comment.char="")
colnames(d.gene) <- c("chr","start","end","attChr","attStart","attEnd","id")
d.gene$feature<-"Gene"

d.ltr<-read.table("bt-intersect.CorFlo.w100000_s50000-LTR.tsv",
                  sep="\t",
                  header=FALSE,
                  check.names=FALSE,
                  comment.char="")
colnames(d.ltr) <- c("chr","start","end","attChr","attStart","attEnd","id")
d.ltr$feature<-"LTR"

d.dna<-read.table("bt-intersect.CorFlo.w100000_s50000-DNA.tsv",
                  sep="\t",
                  header=FALSE,
                  check.names=FALSE,
                  comment.char="")
colnames(d.dna) <- c("chr","start","end","attChr","attStart","attEnd","id")
d.dna$feature<-"DNA transposon"

d.all<-rbind(d.gene, d.ltr, d.dna)

d.sub<-select(d.all, chr, start, end, feature)

#count number of features per window and stardize by max feature per all windows
#this allows us to plot everything on the same y-axis
d.count<-d.sub %>% group_by(feature) %>% count(chr, start) %>% mutate(stdn = n/max(n, na.rm=TRUE))

#rearrange dataframe so contigs are in descending order
d.count$chr<-factor(d.count$chr, levels=c("Contig0", "Contig1", "Contig2", "Contig3",
                                          "Contig4", "Contig5", "Contig6", "Contig7",
                                          "Contig8", "Contig9", "Contig10"))

#divide start values by 1,000,000
ggplot(d.count, aes(x=start/1000000, y=stdn, color=feature)) +
  #use loess smoothing method (to generalize, a moving average)
  #span controls the amount of smoothing (smaller == wigglier lines, larger == smoother lines)
  geom_smooth(span=0.1, method="loess") +
  scale_color_manual(values=c("#0ab45a", "#aa0a3c","#005ac8"),
                     name="Feature",
                     #change order of legend
                     breaks=c("Gene","DNA transposon","LTR")) +
  facet_grid(chr~.) +
  #change x-axis to range from 0-160 Mbp with breaks every 10 Mbp
  scale_x_continuous(breaks=seq(0, 160, 10), name="Chromosome position (Mbp)") +
  scale_y_continuous(name="Relative density") +
  theme_classic()

ggsave("featRelDensity.pdf", height=11, width=8.5, useDingbats=FALSE)
