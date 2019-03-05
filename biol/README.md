# Analysis of biological data

## Differential gene expression

The purpose of this lab is to get a better understanding of how to use the [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) package. I will be using [Bewick, Ji, Niederhuth, Willing et al. (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27457936) RNA-seq data for _Arabidopsis thaliana_ wild type and eighth-generation _met1_ epigenetic Recombinant Inbred Lines (epiRILs) ([Reinders et al. 2009](https://www.ncbi.nlm.nih.gov/pubmed/19390088)). BAM and GTF files are located in `files/`. I have chose to annotate the `R` script `edger.r` in `scripts/` because differential gene expression analysis is fairly involved.

__Putting the data into the right format for edgeR.__ `edgeR` works on a table of integer read counts, with rows corresponding to genes and columns to independent libraries. `edgeR` stores data in a simple list-based data object called a `DGEList`. This type of object is easy to use because it can be manipulated like any list in `R`. You can make this in `R` by specifying the counts and the groups in the function `DGEList()`.

```
setwd("/YOUR/WORKING/DIRECTORY")

library(Rsubread); library(edgeR)

#vector of bam files in current directory
files<-dir(".", "bam")
#quantify read counts for each gene
fc<-featureCounts(files,
                  annot.ext="Athaliana_447_Araport11.gene_exons.gtf",
                  isGTFAnnotationFile=TRUE,
                  GTF.featureType="exon",
                  GTF.attrType="gene_id",
                  isPairedEnd=FALSE,
                  requireBothEndsMapped=TRUE,
                  useMetaFeatures=TRUE)
#vector of files is in alphabetical order
dataGroups<-c("met1", "met1", "met1", "wt", "wt", "wt")
d<-DGEList(counts=cbind(fc$counts), genes=fc$annotation, group=factor(dataGroups))
```

__Filtering teh data.__ First get rid of genes which did not occur frequently enough. Here the cutoff of 0.5 for the CPM has been chosen because it is roughly equal to 10/_L_ where _L_ is the minimum library size in millions. The library sizes here are 18–27 million. We used a round value of 0.5 just for simplicity; the exact value is not important because the downstream differential expression analysis is not sensitive to the small changes in this parameter. The requirement of ≥4 libraries is because each group contains three replicates. This ensures that a gene will be retained if it is expressed in both groups.

```
dim(d)
#[1] 27655     6
#keep original...just in case
d.full<-d
keep<-rowSums(cpm(d)>0.5) >= 4
d<-d[keep,]
dim(d)
#[1] 19815     6
#a reduction of 7840 genes
d$samples$lib.size<-colSums(d$counts)
```

__Normalizing the data.__ `edgeR` is concerned with differential expression analysis rather than with the quantification of expression levels. It is concerned with relative changes in expression levels between conditions, but not directly with estimating absolute expression levels.

The `calcNormFactors()` function normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes. The default method for computing these scale factors uses a trimmed mean of M-values (TMM) between each pair of samples. We call the product of the original library size and the scaling factor the effective library size. The effective library size replaces the original library size in all downsteam analyses.

```
d<-calcNormFactors(d)
```

__Data Exploration.__ Before proceeding with the computations for differential expression, it is possible to produce a plot showing the sample relations based on multidimensional scaling. The basic premise is that we make a plot so samples which are similar are near to each other in the plot while samples that are dissimilar are far from each other.

```
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("topright", as.character(unique(d$samples$group)), col=1:3, pch=20)
#export table of RPKM for further exploration
rpkm<-rpkm(d)
write.table(rpkm, file="rpkm.tsv", quote=F)
```

__Estimating the Dispersion.__ The first major step in the analysis of DGE data using the NB model is to estimate the dispersion parameter for each tag, a measure of the degree of inter-library variation for that tag. Estimating the common dispersion gives an idea of overall variability across the genome for this dataset.

In this example, I am renaming the variable to d1 because we can estimate dispersion by assuming everything has the same common dispersion, or we can use a generalized linear model to try to estimate the dispersion. For now, we will just use the naive method of assuming all tags have the same dispersion.

```
d1<-estimateCommonDisp(d, verbose=T)
d1<-estimateTagwiseDisp(d1)
plotBCV(d1)
```

__GLM estimates of dispersion.__ Fitting a model in `edgeR` takes several steps. First, you must fit the common dispersion. Then you need to fit a trended model (if you do not fit a trend, the default is to use the common dispersion as a trend). Then you can fit the tagwise dispersion which is a function of this model.

In addition to the common and tagwise disperson, we can also estimate a generalized linear model (glm) fit using `edgeR`. In the same way that we've been doing above, we will just add these as additional data to the object we've been working with.

```
design.mat<-model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="auto")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)
```

__Differential Expression.__ Once the dispersions are estimated, we can proceed with testing procedures for determining differential expression. The function `exactTest()` conducts tagwise tests using the exact negative binomial test. The test results for the _n_ most significant tags are conveniently displayed by the `topTags()` function. By default, Benjamini and Hochberg's algorithm is used to control the false discovery rate (FDR).

Recall that d1 is the naive method where we only fit a common dispersion.

```
et<-exactTest(d1, pair=c("met1","wt"))
de1<-decideTestsDGE(et, adjust.method="BH", p.value=0.05)
summary(de1)
#       wt-met1
#Down       819
#NotSig   17644
#Up        1352
```

The function `plotSmear` generates a plot of the tagwise log-fold-changes against log-cpm (analogous to an MA-plot for microarray data). DE tags are highlighted on the plot:

```
de1tags12<-rownames(d1)[as.logical(de1)] 
plotSmear(et, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")
```

__GLM testing for differential expression.__ Just as we used a GLM to fit the trend line above, we can also use this in finding the tags that are interesting by using a likelihood ratio test.

```
fit <- glmFit(d2, design.mat)
lrt<-glmLRT(fit, contrast=c(-1,1))
de2 <- decideTestsDGE(lrt, adjust.method="BH", p.value = 0.05)
summary(de2)
#       -1*met1 1*wt
#Down            813
#NotSig        17639
#Up             1363
de2tags12 <- rownames(d2)[as.logical(de2)]
plotSmear(lrt, de.tags=de2tags12)
abline(h = c(-2, 2), col = "blue")
```

## Gene Ontology (GO) term enrichment

Gene Ontology is a major bioinformatics initiative to unify the representation of gene and gene product attributes across all species. It provides a system for hierarchically classifying genes or gene products to terms organized in a graph structure called an ontology. The terms are grouped into three categories: __Molecular Function (MF)__ (describing the molecular activity of a gene), __Biological Process (BP)__ (describing the larger cellular or physiological role carried out by the gene, coordinated with other genes) and __Cellular Component (CC)__ (describing the location in the cell where the gene product executes its function). Each gene can be described (annotated) with multiple terms. Gene Ontology enrichment involves testing for a biased presence of certain GO terms population of genes (Genes of Interest (GOI)) compared to all genes (universe).

There are several programs to test for GO term enrichment. However, my program of choice is [topGO](http://bioconductor.org/packages/release/bioc/html/topGO.html) mostly because it is straightfoward to use and has good user documentation (See script `topGO.r`). I would avoid drawing too large of conclusions from GO term enrichment...
