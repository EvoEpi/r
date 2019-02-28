library(ggplot2); library(plyr); library(dplyr); library(gridExtra)

setwd("/PATH/TO/YOUR/WORKING/DIRECTORY")

d.col0cg<-read.table("Col-0_CG_per_site_methylation.tsv",
                 sep="\t",
                 header=T,
                 check.names=FALSE,
                 comment.char="")

d.col0chg<-read.table("Col-0_CHG_per_site_methylation.tsv",
                  sep="\t",
                  header=T,
                  check.names=FALSE,
                  comment.char="")

d.col0chh<-read.table("Col-0_CHH_per_site_methylation.tsv",
                  sep="\t",
                  header=T,
                  check.names=FALSE,
                  comment.char="")

d.col0<-rbind(d.col0cg, d.col0chg, d.col0chh)
d.col0$geno<-"wt"

d.met1cg<-read.table("met1_CG_per_site_methylation.tsv",
                     sep="\t",
                     header=T,
                     check.names=FALSE,
                     comment.char="")

d.met1chg<-read.table("met1_CHG_per_site_methylation.tsv",
                      sep="\t",
                      header=T,
                      check.names=FALSE,
                      comment.char="")

d.met1chh<-read.table("met1_CHH_per_site_methylation.tsv",
                      sep="\t",
                      header=T,
                      check.names=FALSE,
                      comment.char="")

d.met1<-rbind(d.met1cg, d.met1chg, d.met1chh)
d.met1$geno<-"met1"

d.all<-rbind(d.col0, d.met1)

p1<-ggplot(d.all, aes(x=per_site_wei, fill=geno)) +
  geom_density(alpha=0.5) +
  facet_grid(con~.) +
  scale_fill_manual(values=c("#aa0a3c","#005ac8"),
                    name="Genotype") +
  ylab("Density") +
  xlab("Per-methylated site weighted methylation") +
  theme_classic()

#count number of methylated sites per genotype per sequence context
d.count<-count(d.all, c("geno", "con"))

p2<-ggplot(d.count, aes(y=freq, x=con, fill=geno)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#aa0a3c","#005ac8"),
                    name="Genotype") +
  ylab("Number of methylated sites") +
  xlab("Sequence context") +
  theme_classic()

#merge plots within one grid (and visualize)
grid.arrange(p1, p2, ncol = 1)

#save
g<-arrangeGrob(p1, p2, ncol=1)
ggsave("grid.pdf", g, height=8.5, width=11, useDingbats=FALSE)
