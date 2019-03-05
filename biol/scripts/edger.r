setwd("/YOUR/WORKING/DIRECTORY")

library(Rsubread); library(edgeR)

files<-dir(".", "bam")
fc<-featureCounts(files,
                  annot.ext="Athaliana_447_Araport11.gene_exons.gtf",
                  isGTFAnnotationFile=TRUE,
                  GTF.featureType="exon",
                  GTF.attrType="gene_id",
                  isPairedEnd=FALSE,
                  requireBothEndsMapped=TRUE,
                  useMetaFeatures=TRUE)
dataGroups<-c("met1", "met1", "met1", "wt", "wt", "wt")
d<-DGEList(counts=cbind(fc$counts), genes=fc$annotation, group=factor(dataGroups))
dim(d)
d.full<-d
keep<-rowSums(cpm(d)>0.5) >= 4
d<-d[keep,]
dim(d)
d$samples$lib.size<-colSums(d$counts)
d<-calcNormFactors(d)
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("topright", as.character(unique(d$samples$group)), col=1:3, pch=20)
rpkm<-rpkm(d)
write.table(rpkm, file="rpkm.tsv", quote=F)
d1<-estimateCommonDisp(d, verbose=T)
d1<-estimateTagwiseDisp(d1)
plotBCV(d1)
design.mat<-model.matrix(~ 0 + d$samples$group)
colnames(design.mat)<-levels(d$samples$group)
d2<-estimateGLMCommonDisp(d,design.mat)
d2<-estimateGLMTrendedDisp(d2,design.mat, method="auto")
d2<-estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)
et<-exactTest(d1, pair=c("met1","wt"))
de1<-decideTestsDGE(et, adjust.method="BH", p.value=0.05)
summary(de1)
de1tags12<-rownames(d1)[as.logical(de1)] 
plotSmear(et, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
fit<-glmFit(d2, design.mat)
lrt<-glmLRT(fit, contrast=c(-1,1))
de2<-decideTestsDGE(lrt, adjust.method="BH", p.value=0.05)
summary(de2)
de2tags<-rownames(d2)[as.logical(de2)]
plotSmear(lrt, de.tags=de2tags)
abline(h = c(-2, 2), col="blue")
