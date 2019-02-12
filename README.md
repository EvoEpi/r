# r
A collection of R scripts for analysis of biological data, plotting, and stats.

## Analyses

__Gene Ontology (GO) term enrichment__. Gene Ontology is a major bioinformatics initiative to unify the representation of gene and gene product attributes across all species. It provides a system for hierarchically classifying genes or gene products to terms organized in a graph structure called an ontology. The terms are grouped into three categories: __Molecular Function (MF)__ (describing the molecular activity of a gene), __Biological Process (BP)__ (describing the larger cellular or physiological role carried out by the gene, coordinated with other genes) and __Cellular Component (CC)__ (describing the location in the cell where the gene product executes its function). Each gene can be described (annotated) with multiple terms. Gene Ontology enrichment involves testing for a biased presence of certain GO terms population of genes (Genes of Interest (GOI)) compared to all genes (universe).

There are several programs to test for GO term enrichment. However, my program of choice is [topGO](http://bioconductor.org/packages/release/bioc/html/topGO.html) mostly because it is straightfoward to use and has good user documentation (See scrip `topGO.r`). I would avoid drawing too large of conclusions from GO term enrichment...

## Stats

__Permutation tests (as described by [wikipedia](https://en.wikipedia.org/wiki/Resampling_(statistics)#Permutation_tests))__. A permutation test is a type of statistical significance test in which the distribution of the test statistic under the null hypothesis is obtained by calculating all possible values of the test statistic under rearrangements of the labels on the observed data points. In other words, the method by which treatments are allocated to subjects in an experimental design is mirrored in the analysis of that design. If the labels are exchangeable under the null hypothesis, then the resulting tests yield exact significance levels. Confidence intervals can then be derived from the tests. The theory has evolved from the works of Ronald Fisher and EJG Pitman in the 1930s.

I have generated two scripts – `permutationTestDist.r` and `permutationTestMean.r` – as laid out by this [StackExchange](https://stats.stackexchange.com/questions/136661/using-bootstrap-under-h0-to-perform-a-test-for-the-difference-of-two-means-repl) discussion to test whether the the distribution of _x_ and _y_ are identical and whether or not their population means are equal, without making any assumptions about their variance. To use, load your list of values in _x_ and _y_.
