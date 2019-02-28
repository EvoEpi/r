library(ggplot2)

setwd("/PATH/TO/YOUR/WORKING/DIRECTORY")

d.col0<-read.table("Col-0_metaplot.tsv",
                 sep="\t",
                 header=T,
                 check.names=FALSE,
                 comment.char="")
d.col0$geno<-"wt"

d.met1<-read.table("met1_metaplot.tsv",
                  sep="\t",
                  header=T,
                  check.names=FALSE,
                  comment.char="")
d.met1$geno<-"met1"

d.all<-rbind(d.col0, d.met1)

ggplot(d.all, aes(y=wei, x=bin, color=geno)) +
  geom_line() +
  facet_grid(con~., scales="free") +
  scale_color_manual(values=c("#aa0a3c","#005ac8"),
                    name="Genotype") +
  ylab("Weighted methylation") +
  xlab("Bin") +
  scale_x_continuous(limits=c(0, 60),
                    #limits and breaks should correspond
                    #to input values of metaplot.py
                    breaks=c(0, 20, 40, 60),
                    #flanking region should correspond
                    #to input values of metaplot.py
                    labels=c("-1kb", "TSS", "TTS", "+1kb"),
                    name="Position") +
  theme_classic()
